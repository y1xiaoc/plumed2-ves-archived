/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2011-2015 The plumed team
   (see the PEOPLE file at the root of the distribution for a list of names)

   See http://www.plumed-code.org for more information.

   This file is part of plumed, version 2.

   plumed is free software: you can redistribute it and/or modify
   it under the terms of the GNU Lesser General Public License as published by
   the Free Software Foundation, either version 3 of the License, or
   (at your option) any later version.

   plumed is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU Lesser General Public License for more details.

   You should have received a copy of the GNU Lesser General Public License
   along with plumed.  If not, see <http://www.gnu.org/licenses/>.
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
#include "VesBias.h"
#include "ves_basisfunctions/BasisFunctions.h"
#include "ves_tools/CoeffsVector.h"
#include "ves_tools/CoeffsMatrix.h"

#include "tools/Communicator.h"
#include "core/PlumedMain.h"
#include "core/Atoms.h"


namespace PLMD{
namespace bias{

VesBias::VesBias(const ActionOptions&ao):
Bias(ao),
coeffs_ptr(NULL),
coeffderivs_aver_ps_ptr(NULL),
gradient_ptr(NULL),
hessian_ptr(NULL),
coeffderivs_aver_sampled(0),
hessian_diagonal_(true),
use_mwalkers_(false),
aver_counter(0.0),
kbt_(0.0),
beta_(0.0)
{
  bool full_hessian=false;
  parseFlag("FULL_HESSIAN",full_hessian);
  hessian_diagonal_ = !full_hessian;
  //
  parseFlag("MULTIPLE_WALKERS",use_mwalkers_);
  std::cerr << use_mwalkers_;
  if(use_mwalkers_){
   log.printf("  Using multiple walkers:\n");
   log.printf("   number of walkers: %d\n",multi_sim_comm.Get_size());
   log.printf("   walker number: %d\n",multi_sim_comm.Get_rank());
 }
 double temp=0.0;
 parse("TEMP",temp);
 if(temp>0.0) kbt_=plumed.getAtoms().getKBoltzmann()*temp;
 else kbt_=plumed.getAtoms().getKbT();
 if(kbt_==0.0) error("the MD engine does not pass the temperature to plumed so you need to give with the TEMP keyword");
 beta_=1.0/kbt_;
}

VesBias::~VesBias(){
  delete coeffs_ptr;
  delete coeffderivs_aver_ps_ptr;
  delete gradient_ptr;
  delete hessian_ptr;
}


void VesBias::registerKeywords( Keywords& keys ) {
  Bias::registerKeywords(keys);
  keys.addFlag("FULL_HESSIAN",false,"if the full Hessian matrix should be used for the optimization, otherwise only the diagonal Hessian is used");
  keys.addFlag("MULTIPLE_WALKERS",false,"if optimization is to be performed using multiple walkers connected via MPI");
  keys.add("optional","TEMP","the system temperature - this is needed if the MD code does not pass the temperature");
}


void VesBias::initializeCoeffs(const std::vector<std::string>& dimension_labels,const std::vector<unsigned int>& indices_shape) {
  coeffs_ptr = new CoeffsVector("coeffs",dimension_labels,indices_shape,comm,true);
  initializeGradientAndHessian();
}


void VesBias::initializeCoeffs(std::vector<Value*>& args,std::vector<BasisFunctions*>& basisf) {
  coeffs_ptr = new CoeffsVector("coeffs",args,basisf,comm,true);
  initializeGradientAndHessian();
}


void VesBias::linkCoeffs(CoeffsVector* coeffs_ptr_in) {
  coeffs_ptr = coeffs_ptr_in;
  initializeGradientAndHessian();
}


void VesBias::linkCoeffs(CoeffsVector& coeffs_in) {
  coeffs_ptr = &coeffs_in;
  initializeGradientAndHessian();
}


void VesBias::initializeGradientAndHessian() {
  //
  coeffderivs_aver_ps_ptr = new CoeffsVector(*coeffs_ptr);
  coeffderivs_aver_ps_ptr->setLabels("average-over-ps");
  //
  gradient_ptr = new CoeffsVector(*coeffs_ptr);
  gradient_ptr->setLabels("gradient");
  //
  hessian_ptr = new CoeffsMatrix("hessian",coeffs_ptr,comm,hessian_diagonal_,false);
  //
  coeffderivs_cov_sampled.assign(hessian_ptr->getSize(),0.0);
  coeffderivs_aver_sampled.assign(numberOfCoeffs(),0.0);
  aver_counter=0.0;
}


void VesBias::updateGradientAndHessian() {
  comm.Sum(coeffderivs_aver_sampled);
  comm.Sum(coeffderivs_cov_sampled);
  Gradient() = CoeffDerivsAverTargetDist() - coeffderivs_aver_sampled;
  Hessian() = coeffderivs_cov_sampled;
  Hessian() *= getBeta();
  if(use_mwalkers_){
    Gradient().sumMultiSimCommMPI(multi_sim_comm);
    Hessian().sumMultiSimCommMPI(multi_sim_comm);
  }
  //
  coeffderivs_cov_sampled.assign(hessian_ptr->getSize(),0.0);
  coeffderivs_aver_sampled.assign(numberOfCoeffs(),0.0);
  aver_counter=0.0;
}


void VesBias::clearGradientAndHessian() {}


void VesBias::setCoeffsDerivs(const std::vector<double>& coeffderivs) {
  /*
  use the following online equation to calculate the average and covariance (see wikipedia)
      xm[n+1] = xm[n] + (x[n+1]-xm[n])/(n+1)
      cov(x,y)[n+1] = ( cov(x,y)[n]*n + (n/(n+1))*(x[n+1]-xm[n])*(y[n+1]-ym[n]) ) / (n+1)
                    = cov(x,y)[n]*(n/(n+1)) + ( n * (x[n+1]-xm[n])/(n+1) * (y[n+1]-ym[n])/(n+1) );
      n starts at 0.
  */
  size_t ncoeffs = numberOfCoeffs();
  std::vector<double> deltas(ncoeffs,0.0);
  size_t stride = comm.Get_size();
  size_t rank = comm.Get_rank();
  // update average and diagonal part of Hessian
  for(size_t i=rank; i<ncoeffs;i+=stride){
    size_t midx = hessian_ptr->getMatrixIndex(i,i);
    deltas[i] = (coeffderivs[i]-coeffderivs_aver_sampled[i])/(aver_counter+1); // (x[n+1]-xm[n])/(n+1)
    coeffderivs_aver_sampled[i] += deltas[i];
    coeffderivs_cov_sampled[midx] = coeffderivs_cov_sampled[midx] * ( aver_counter / (aver_counter+1) ) + aver_counter*deltas[i]*deltas[i];
  }
  comm.Sum(deltas);
  // update off diagonal part of Hessian
  if(!hessian_diagonal_){
    for(size_t i=rank; i<ncoeffs;i+=stride){
      for(size_t j=(i+1); j<ncoeffs;j++){
        size_t midx = hessian_ptr->getMatrixIndex(i,j);
        coeffderivs_cov_sampled[midx] = coeffderivs_cov_sampled[midx] * ( aver_counter / (aver_counter+1) ) + aver_counter*deltas[i]*deltas[j];
      }
    }
  }
  // NOTE: the MPI sum for coeffderivs_aver_sampled and coeffderivs_cov_sampled is done later
  aver_counter += 1.0;
}

void VesBias::setCoeffsDerivsOverTargetDist(const std::vector<double>& coeffderivs_aver_ps) {
  CoeffDerivsAverTargetDist() = coeffderivs_aver_ps;
}


}
}
