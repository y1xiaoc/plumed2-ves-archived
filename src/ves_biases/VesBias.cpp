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
#include "ves_optimizers/Optimizer.h"

#include "tools/Communicator.h"
#include "core/PlumedMain.h"
#include "core/Atoms.h"

#include <iostream>


namespace PLMD{
namespace bias{

VesBias::VesBias(const ActionOptions&ao):
Bias(ao),
ncoeffssets_(0),
coeffs_pntrs(0),
coeffderivs_aver_ps_pntrs(0),
gradient_pntrs(0),
hessian_pntrs(0),
coeffderivs_aver_sampled(0),
coeffderivs_cov_sampled(0),
use_multiple_coeffssets_(false),
ncoeffs_total_(0),
optimizer_pntr(NULL),
optimize_coeffs_(false),
compute_hessian_(false),
diagonal_hessian_(true),
aver_counter(0.0),
kbt_(0.0)
{
  double temp=0.0;
  parse("TEMP",temp);
  if(temp>0.0){
    kbt_=plumed.getAtoms().getKBoltzmann()*temp;
  }
  else {
    kbt_=plumed.getAtoms().getKbT();
  }
  if(kbt_>0.0){
    log.printf("  KbT: %f\n",kbt_);
  }
  // NOTE: the check for that the temperature is given is done when linking the optimizer later on.
}

VesBias::~VesBias(){
  for(unsigned int i=0; i<ncoeffssets_; i++){
    delete coeffs_pntrs[i];
    delete coeffderivs_aver_ps_pntrs[i];
    delete gradient_pntrs[i];
    delete hessian_pntrs[i];
  }
}


void VesBias::registerKeywords( Keywords& keys ) {
  Bias::registerKeywords(keys);
  keys.add("optional","TEMP","the system temperature - this is needed if the MD code does not pass the temperature to PLUMED");
  keys.addOutputComponent("bias","default","the instantaneous value of the bias potential");
  keys.addOutputComponent("force2","default","the instantaneous value of the squared force due to this bias potential");
}


void VesBias::addCoeffsSet(const std::vector<std::string>& dimension_labels,const std::vector<unsigned int>& indices_shape) {
  CoeffsVector* coeffs_pntr_tmp = new CoeffsVector("coeffs",dimension_labels,indices_shape,comm,true);
  coeffs_pntr_tmp->linkVesBias(this);
  initializeCoeffs(coeffs_pntr_tmp);
}


void VesBias::addCoeffsSet(std::vector<Value*>& args,std::vector<BasisFunctions*>& basisf) {
  CoeffsVector* coeffs_pntr_tmp = new CoeffsVector("coeffs",args,basisf,comm,true);
  coeffs_pntr_tmp->linkVesBias(this);
  initializeCoeffs(coeffs_pntr_tmp);
}


void VesBias::initializeCoeffs(CoeffsVector* coeffs_pntr_in) {
  //
  std::string label;
  if(!use_multiple_coeffssets_ && ncoeffssets_==1){
    plumed_merror("you are not allowed to use multiple coefficient sets");
  }
  //
  label = labelString("coeffs",ncoeffssets_);
  coeffs_pntr_in->setLabels(label);

  coeffs_pntrs.push_back(coeffs_pntr_in);
  CoeffsVector* aver_ps_tmp = new CoeffsVector(*coeffs_pntr_in);
  label = labelString("ps-aver",ncoeffssets_);
  aver_ps_tmp->setLabels(label);
  aver_ps_tmp->setValues(0.0);
  coeffderivs_aver_ps_pntrs.push_back(aver_ps_tmp);
  //
  CoeffsVector* gradient_tmp = new CoeffsVector(*coeffs_pntr_in);
  label = labelString("gradient",ncoeffssets_);
  gradient_tmp->setLabels(label);
  gradient_pntrs.push_back(gradient_tmp);
  //
  label = labelString("hessian",ncoeffssets_);
  CoeffsMatrix* hessian_tmp = new CoeffsMatrix(label,coeffs_pntr_in,comm,diagonal_hessian_,true);
  hessian_pntrs.push_back(hessian_tmp);
  //
  std::vector<double> aver_sampled_tmp;
  aver_sampled_tmp.assign(coeffs_pntr_in->numberOfCoeffs(),0.0);
  coeffderivs_aver_sampled.push_back(aver_sampled_tmp);
  //
  std::vector<double> cov_sampled_tmp;
  cov_sampled_tmp.assign(hessian_tmp->getSize(),0.0);
  coeffderivs_cov_sampled.push_back(cov_sampled_tmp);
  //
  ncoeffssets_++;
}


void VesBias::updateGradientAndHessian() {
  for(unsigned int i=0; i<ncoeffssets_; i++){
    comm.Sum(coeffderivs_aver_sampled[i]);
    comm.Sum(coeffderivs_cov_sampled[i]);
    Gradient(i) = CoeffDerivsAverTargetDist(i) - coeffderivs_aver_sampled[i];
    Hessian(i) = coeffderivs_cov_sampled[i];
    Hessian(i) *= getBeta();
    std::fill(coeffderivs_aver_sampled[i].begin(), coeffderivs_aver_sampled[i].end(), 0.0);
    std::fill(coeffderivs_cov_sampled[i].begin(), coeffderivs_cov_sampled[i].end(), 0.0);
  }
  aver_counter=0.0;
}


void VesBias::clearGradientAndHessian() {}


void VesBias::setCoeffsDerivs(const std::vector<double>& coeffderivs, const unsigned c_id) {
  /*
  use the following online equation to calculate the average and covariance (see wikipedia)
      xm[n+1] = xm[n] + (x[n+1]-xm[n])/(n+1)
      cov(x,y)[n+1] = ( cov(x,y)[n]*n + (n/(n+1))*(x[n+1]-xm[n])*(y[n+1]-ym[n]) ) / (n+1)
                    = cov(x,y)[n]*(n/(n+1)) + ( n * (x[n+1]-xm[n])/(n+1) * (y[n+1]-ym[n])/(n+1) );
      n starts at 0.
  */
  size_t ncoeffs = numberOfCoeffs(c_id);
  std::vector<double> deltas(ncoeffs,0.0);
  size_t stride = comm.Get_size();
  size_t rank = comm.Get_rank();
  // update average and diagonal part of Hessian
  for(size_t i=rank; i<ncoeffs;i+=stride){
    size_t midx = getHessianIndex(i,i,c_id);
    deltas[i] = (coeffderivs[i]-coeffderivs_aver_sampled[c_id][i])/(aver_counter+1); // (x[n+1]-xm[n])/(n+1)
    coeffderivs_aver_sampled[c_id][i] += deltas[i];
    coeffderivs_cov_sampled[c_id][midx] = coeffderivs_cov_sampled[c_id][midx] * ( aver_counter / (aver_counter+1) ) + aver_counter*deltas[i]*deltas[i];
  }
  comm.Sum(deltas);
  // update off-diagonal part of the Hessian
  if(!diagonal_hessian_){
    for(size_t i=rank; i<ncoeffs;i+=stride){
      for(size_t j=(i+1); j<ncoeffs;j++){
        size_t midx = getHessianIndex(i,j,c_id);
        coeffderivs_cov_sampled[c_id][midx] = coeffderivs_cov_sampled[c_id][midx] * ( aver_counter / (aver_counter+1) ) + aver_counter*deltas[i]*deltas[j];
      }
    }
  }
  // NOTE: the MPI sum for coeffderivs_aver_sampled and coeffderivs_cov_sampled is done later
  //aver_counter += 1.0;
}


void VesBias::setCoeffsDerivsOverTargetDist(const std::vector<double>& coeffderivs_aver_ps, const unsigned coeffs_id) {
  CoeffDerivsAverTargetDist(coeffs_id) = coeffderivs_aver_ps;
}


void VesBias::linkOptimizer(Optimizer* optimizer_pntr_in) {
  //
  if(optimizer_pntr==NULL){
    optimizer_pntr = optimizer_pntr_in;
  }
  else {
    std::string err_msg = "VES bias " + getName() + " with label " + getLabel() + " has already been linked with optimizer " + optimizer_pntr->getName() + " with label " + optimizer_pntr->getLabel() + ". You cannot link two optimizer to the same VES bias.";
    plumed_merror(err_msg);
  }
  //
  if(kbt_==0.0){
    std::string err_msg = "VES bias " + getName() + " with label " + getLabel() + ": if you want to optimize this bias you need to give the temperature using the TEMP keyword as the MD engine does not pass it to PLUMED";
    plumed_merror(err_msg);
  }
  //
  optimize_coeffs_ = true;
}


void VesBias::enableHessian(const bool diagonal_hessian) {
  compute_hessian_=true;
  diagonal_hessian_=diagonal_hessian;
  coeffderivs_cov_sampled.clear();
  for (unsigned int i=0; i<ncoeffssets_; i++){
    delete hessian_pntrs[i];
    std::string label = labelString("hessian",i);
    hessian_pntrs[i] = new CoeffsMatrix(label,coeffs_pntrs[i],comm,diagonal_hessian_,true);
    std::vector<double> cov_sampled_tmp;
    cov_sampled_tmp.assign(hessian_pntrs[i]->getSize(),0.0);
    coeffderivs_cov_sampled.push_back(cov_sampled_tmp);
  }
}


void VesBias::disableHessian() {
  compute_hessian_=false;
  diagonal_hessian_=true;
  coeffderivs_cov_sampled.clear();
  for (unsigned int i=0; i<ncoeffssets_; i++){
    delete hessian_pntrs[i];
    std::string label = labelString("hessian",i);
    hessian_pntrs[i] = new CoeffsMatrix(label,coeffs_pntrs[i],comm,diagonal_hessian_,true);
    std::vector<double> cov_sampled_tmp;
    cov_sampled_tmp.assign(hessian_pntrs[i]->getSize(),0.0);
    coeffderivs_cov_sampled.push_back(cov_sampled_tmp);
  }
}


void VesBias::apply() {
  Bias::apply();
  aver_counter += 1.0;
}


std::string VesBias::labelString(const std::string& type, const unsigned int coeffs_id) {
  std::string label_prefix = getLabel() + ".";
  std::string label_postfix = "";
  if(use_multiple_coeffssets_){
    Tools::convert(coeffs_id,label_postfix);
    label_postfix = "-" + label_postfix;
  }
  return label_prefix+type+label_postfix;
}


}
}
