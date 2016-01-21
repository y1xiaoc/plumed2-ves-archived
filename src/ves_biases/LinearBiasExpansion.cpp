/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2012-2014 The plumed team
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
#include "LinearBiasExpansion.h"
#include "VesBias.h"
#include "ves_tools/CoeffsVector.h"
#include "ves_basisfunctions/BasisFunctions.h"


#include "tools/Keywords.h"
#include "tools/Grid.h"
#include "tools/Communicator.h"


namespace PLMD{
namespace bias{

//+PLUMEDOC VES LinearBiasExpansion
/*
*/
//+ENDPLUMEDOC

void LinearBiasExpansion::registerKeywords(Keywords& keys){
}


LinearBiasExpansion::LinearBiasExpansion(
  const std::string& label,
  Communicator& cc,
  std::vector<Value*> args_pntrs_in,
  std::vector<BasisFunctions*> basisf_pntrs_in,
  CoeffsVector* bias_coeffs_pntr_in):
label_(label),
action_pntr(NULL),
vesbias_pntr(NULL),
mycomm(cc),
serial_(false),
args_pntrs(args_pntrs_in),
nargs_(args_pntrs.size()),
basisf_pntrs(basisf_pntrs_in),
nbasisf_(nargs_),
bias_coeffs_pntr(NULL),
ncoeffs_(0),
coeffderivs_aver_ps_pntr(NULL),
fes_wt_coeffs_pntr(NULL),
biasf_(-1.0),
invbiasf_(-1.0),
bias_grid_pntr(NULL),
fes_grid_pntr(NULL),
ps_grid_pntr(NULL)
{
  plumed_massert(args_pntrs.size()==basisf_pntrs.size(),"number of arguments and basis functions do not match");
  for(unsigned int k=0;k<nargs_;k++){nbasisf_[k]=basisf_pntrs[k]->getNumberOfBasisFunctions();}
  //
  if(bias_coeffs_pntr_in==NULL){
    bias_coeffs_pntr = new CoeffsVector(label_+".coeffs",args_pntrs,basisf_pntrs,mycomm,true);
  }
  //
  coeffderivs_aver_ps_pntr = new CoeffsVector(*bias_coeffs_pntr);
  std::string coeffderivs_aver_ps_label = bias_coeffs_pntr->getLabel();
  if(coeffderivs_aver_ps_label.find("coeffs")!=std::string::npos){
    coeffderivs_aver_ps_label.replace(coeffderivs_aver_ps_label.find("coeffs"), std::string("coeffs").length(), "coeffderivs_aver_ps");
  }
  else {
    coeffderivs_aver_ps_label += "_aver_ps";
  }
  coeffderivs_aver_ps_pntr->setLabels(coeffderivs_aver_ps_label);
  //
}

LinearBiasExpansion::~LinearBiasExpansion() {
  if(bias_grid_pntr!=NULL){
    delete bias_grid_pntr;
  }
  if(fes_grid_pntr!=NULL){
    delete fes_grid_pntr;
  }
  if(ps_grid_pntr!=NULL){
    delete ps_grid_pntr;
  }
  if(coeffderivs_aver_ps_pntr!=NULL){
    delete coeffderivs_aver_ps_pntr;
  }
  if(fes_wt_coeffs_pntr!=NULL){
    delete fes_wt_coeffs_pntr;
  }
}


void LinearBiasExpansion::linkVesBias(bias::VesBias* vesbias_pntr_in){
  vesbias_pntr = vesbias_pntr_in;
  action_pntr = static_cast<Action*>(vesbias_pntr_in);
}


void LinearBiasExpansion::linkAction(Action* action_pntr_in){
  action_pntr = action_pntr_in;
}


void LinearBiasExpansion::setupGrid(const std::vector<unsigned int>& nbins, const bool usederiv){
  plumed_assert(nbins.size()==nargs_);
  std::vector<std::string> min(nargs_);
  std::vector<std::string> max(nargs_);
  for(unsigned int k=0;k<nargs_;k++){
    Tools::convert(basisf_pntrs[k]->intervalMin(),min[k]);
    Tools::convert(basisf_pntrs[k]->intervalMax(),max[k]);
  }
  bias_grid_pntr = new Grid(label_+".bias",args_pntrs,min,max,nbins,false,usederiv);
}


void LinearBiasExpansion::updateBiasGrid(){
  for(unsigned int l=0; l<bias_grid_pntr->getSize(); l++){
    std::vector<double> forces(nargs_);
    std::vector<double> cv_value(nargs_);
    cv_value=bias_grid_pntr->getPoint(l);
    double bias_value=getBiasAndForces(cv_value,forces);
    if(bias_grid_pntr->hasDerivatives()){
      bias_grid_pntr->setValueAndDerivatives(l,bias_value,forces);
    }
    else{
      bias_grid_pntr->setValue(l,bias_value);
    }

 }
}


void LinearBiasExpansion::writeBiasGridToFile(const std::string& filepath, const bool append_file){
  OFile file;
  if(append_file){file.enforceRestart();}
  if(action_pntr!=NULL){
    file.link(*action_pntr);
  }
  file.open(filepath);
  bias_grid_pntr->writeToFile(file);
  file.close();
}


double LinearBiasExpansion::getBiasAndForces(const std::vector<double>& cv_values, std::vector<double>& forces){
  std::vector<double> cv_values_trsfrm(nargs_);
  std::vector<bool>   inside_interval(nargs_,true);
  //
  std::vector< std::vector <double> > bf_values;
  std::vector< std::vector <double> > bf_derivs;
  //
  for(unsigned int k=0;k<nargs_;k++){
    std::vector<double> tmp_val(nbasisf_[k]);
    std::vector<double> tmp_der(nbasisf_[k]);
    bool inside=true;
    basisf_pntrs[k]->getAllValues(cv_values[k],cv_values_trsfrm[k],inside,tmp_val,tmp_der);
    inside_interval[k]=inside;
    bf_values.push_back(tmp_val);
    bf_derivs.push_back(tmp_der);
    forces[k]=0.0;
  }
  //
  unsigned int stride=1;
  unsigned int rank=0;
  if(!serial_)
  {
    stride=mycomm.Get_size();
    rank=mycomm.Get_rank();
  }
  // loop over coeffs
  double bias=0.0;
  for(unsigned int i=rank;i<bias_coeffs_pntr->getSize();i+=stride){
    std::vector<unsigned int> indices=bias_coeffs_pntr->getIndices(i);
    double coeff = bias_coeffs_pntr->getValue(i);
    double bf_curr=1.0;
    for(unsigned int k=0;k<nargs_;k++){bf_curr*=bf_values[k][indices[k]];}
    bias+=coeff*bf_curr;
    for(unsigned int k=0;k<nargs_;k++){
      forces[k]-=coeff*bf_curr*(bf_derivs[k][indices[k]]/bf_values[k][indices[k]]);
    }
  }
  //
  if(!serial_){
    mycomm.Sum(bias);
    mycomm.Sum(forces);
  }
  return bias;
}


double LinearBiasExpansion::getBias(const std::vector<double>& cv_values) {
  std::vector<double> forces(nargs_);
  return getBiasAndForces(cv_values,forces);
}


void LinearBiasExpansion::setupWellTempered(const double biasf, const std::vector<unsigned int>& nbins) {
  plumed_massert(biasf>1.0,"the value of the bias factor doesn't make sense, it should be larger than 1.0");
  biasf_=biasf;
  invbiasf_ = 1.0/biasf_;
  fes_wt_coeffs_pntr = new CoeffsVector(*bias_coeffs_pntr);
  std::string fes_wt_label = bias_coeffs_pntr->getLabel();
  if(fes_wt_label.find("coeffs")!=std::string::npos){
    fes_wt_label.replace(fes_wt_label.find("coeffs"), std::string("coeffs").length(), "fes_wt_coeffs");
  }
  else {
    fes_wt_label += "_fes_wt";
  }
  fes_wt_coeffs_pntr->setLabels(fes_wt_label);
}


void LinearBiasExpansion::updateWellTemperedFESCoeffs() {
  plumed_assert(biasf_>1.0);
  FesWTCoeffs() = -BiasCoeffs() + invbiasf_*FesWTCoeffs();
}


}
}
