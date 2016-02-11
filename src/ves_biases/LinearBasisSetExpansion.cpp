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
#include "LinearBasisSetExpansion.h"
#include "VesBias.h"
#include "ves_tools/CoeffsVector.h"
#include "ves_basisfunctions/BasisFunctions.h"
#include "ves_targetdistributions/TargetDistribution.h"
#include "ves_targetdistributions/TargetDistributionRegister.h"


#include "tools/Keywords.h"
#include "tools/Grid.h"
#include "tools/Communicator.h"

// #include <iostream>



namespace PLMD{
namespace bias{

//+PLUMEDOC VES LinearBasisSetExpansion
/*
*/
//+ENDPLUMEDOC

void LinearBasisSetExpansion::registerKeywords(Keywords& keys) {
}


LinearBasisSetExpansion::LinearBasisSetExpansion(
  const std::string& label,
  Communicator& cc,
  std::vector<Value*> args_pntrs_in,
  std::vector<BasisFunctions*> basisf_pntrs_in,
  CoeffsVector* bias_coeffs_pntr_in):
label_(label),
action_pntr_(NULL),
vesbias_pntr_(NULL),
mycomm_(cc),
serial_(false),
args_pntrs_(args_pntrs_in),
nargs_(args_pntrs_.size()),
basisf_pntrs_(basisf_pntrs_in),
nbasisf_(basisf_pntrs_.size()),
bias_coeffs_pntr_(bias_coeffs_pntr_in),
ncoeffs_(0),
coeffderivs_aver_ps_pntr_(NULL),
fes_wt_coeffs_pntr_(NULL),
biasf_(-1.0),
invbiasf_(-1.0),
bias_grid_pntr_(NULL),
fes_grid_pntr_(NULL),
ps_grid_pntr_(NULL)
{
  plumed_massert(args_pntrs_.size()==basisf_pntrs_.size(),"number of arguments and basis functions do not match");
  for(unsigned int k=0;k<nargs_;k++){nbasisf_[k]=basisf_pntrs_[k]->getNumberOfBasisFunctions();}
  //
  if(bias_coeffs_pntr_==NULL){
    bias_coeffs_pntr_ = new CoeffsVector(label_+".coeffs",args_pntrs_,basisf_pntrs_,mycomm_,true);
  }
  plumed_massert(bias_coeffs_pntr_->numberOfDimensions()==basisf_pntrs_.size(),"dimension of coeffs does not match with number of basis functions ");
  //
  ncoeffs_ = bias_coeffs_pntr_->numberOfCoeffs();
  coeffderivs_aver_ps_pntr_ = new CoeffsVector(*bias_coeffs_pntr_);

  std::string coeffderivs_aver_ps_label = bias_coeffs_pntr_->getLabel();
  if(coeffderivs_aver_ps_label.find("coeffs")!=std::string::npos){
    coeffderivs_aver_ps_label.replace(coeffderivs_aver_ps_label.find("coeffs"), std::string("coeffs").length(), "coeffderivs_aver_ps");
  }
  else {
    coeffderivs_aver_ps_label += "_aver_ps";
  }
  coeffderivs_aver_ps_pntr_->setLabels(coeffderivs_aver_ps_label);
  //
}

LinearBasisSetExpansion::~LinearBasisSetExpansion() {
  if(bias_grid_pntr_!=NULL){
    delete bias_grid_pntr_;
  }
  if(fes_grid_pntr_!=NULL){
    delete fes_grid_pntr_;
  }
  if(ps_grid_pntr_!=NULL){
    delete ps_grid_pntr_;
  }
  if(coeffderivs_aver_ps_pntr_!=NULL){
    delete coeffderivs_aver_ps_pntr_;
  }
  if(fes_wt_coeffs_pntr_!=NULL){
    delete fes_wt_coeffs_pntr_;
  }
}


void LinearBasisSetExpansion::linkVesBias(bias::VesBias* vesbias_pntr_in) {
  vesbias_pntr_ = vesbias_pntr_in;
  action_pntr_ = static_cast<Action*>(vesbias_pntr_in);
}


void LinearBasisSetExpansion::linkAction(Action* action_pntr_in) {
  action_pntr_ = action_pntr_in;
}


void LinearBasisSetExpansion::setupBiasGrid(const std::vector<unsigned int>& nbins, const bool usederiv) {
  plumed_assert(nbins.size()==nargs_);
  std::vector<std::string> min(nargs_);
  std::vector<std::string> max(nargs_);
  for(unsigned int k=0;k<nargs_;k++){
    Tools::convert(basisf_pntrs_[k]->intervalMin(),min[k]);
    Tools::convert(basisf_pntrs_[k]->intervalMax(),max[k]);
  }
  bias_grid_pntr_ = new Grid(label_+".bias",args_pntrs_,min,max,nbins,false,usederiv);
}


void LinearBasisSetExpansion::updateBiasGrid() {
  for(unsigned int l=0; l<bias_grid_pntr_->getSize(); l++){
    std::vector<double> forces(nargs_);
    std::vector<double> coeffsderivs_values(ncoeffs_);
    std::vector<double> args_values = bias_grid_pntr_->getPoint(l);
    double bias_value=getBiasAndForces(args_values,forces,coeffsderivs_values);
    if(bias_grid_pntr_->hasDerivatives()){
      bias_grid_pntr_->setValueAndDerivatives(l,bias_value,forces);
    }
    else{
      bias_grid_pntr_->setValue(l,bias_value);
    }

 }
}


void LinearBasisSetExpansion::writeBiasGridToFile(const std::string& filepath, const bool append_file) {
  OFile file;
  if(append_file){file.enforceRestart();}
  if(action_pntr_!=NULL){
    file.link(*action_pntr_);
  }
  file.open(filepath);
  bias_grid_pntr_->writeToFile(file);
  file.close();
}


double LinearBasisSetExpansion::getBiasAndForces(const std::vector<double>& args_values, std::vector<double>& forces, std::vector<double>& coeffsderivs_values, std::vector<BasisFunctions*>& basisf_pntrs_in, CoeffsVector* coeffs_pntr_in, Communicator* comm_in) {
  unsigned int nargs = args_values.size();
  plumed_assert(coeffs_pntr_in->numberOfDimensions()==nargs);
  plumed_assert(basisf_pntrs_in.size()==nargs);
  plumed_assert(forces.size()==nargs);

  std::vector<double> args_values_trsfrm(nargs);
  std::vector<bool>   inside_interval(nargs,true);
  //
  std::vector< std::vector <double> > bf_values;
  std::vector< std::vector <double> > bf_derivs;
  //
  for(unsigned int k=0;k<nargs;k++){
    std::vector<double> tmp_val(basisf_pntrs_in[k]->getNumberOfBasisFunctions());
    std::vector<double> tmp_der(tmp_val.size());
    bool inside=true;
    basisf_pntrs_in[k]->getAllValues(args_values[k],args_values_trsfrm[k],inside,tmp_val,tmp_der);
    inside_interval[k]=inside;
    bf_values.push_back(tmp_val);
    bf_derivs.push_back(tmp_der);
    forces[k]=0.0;
  }
  //
  size_t stride=1;
  size_t rank=0;
  if(comm_in!=NULL)
  {
    stride=comm_in->Get_size();
    rank=comm_in->Get_rank();
  }
  // loop over coeffs
  double bias=0.0;
  for(size_t i=rank;i<coeffs_pntr_in->numberOfCoeffs();i+=stride){
    std::vector<unsigned int> indices=coeffs_pntr_in->getIndices(i);
    double coeff = coeffs_pntr_in->getValue(i);
    double bf_curr=1.0;
    for(unsigned int k=0;k<nargs;k++){
      bf_curr*=bf_values[k][indices[k]];
    }
    bias+=coeff*bf_curr;
    coeffsderivs_values[i] = bf_curr;
    for(unsigned int k=0;k<nargs;k++){
      forces[k]-=coeff*bf_curr*(bf_derivs[k][indices[k]]/bf_values[k][indices[k]]);
    }
  }
  //
  if(comm_in!=NULL){
    comm_in->Sum(bias);
    comm_in->Sum(forces);
  }
  return bias;
}


void LinearBasisSetExpansion::getBasisSetValues(const std::vector<double>& args_values, std::vector<double>& basisset_values, std::vector<BasisFunctions*>& basisf_pntrs_in, CoeffsVector* coeffs_pntr_in, Communicator* comm_in) {
  unsigned int nargs = args_values.size();
  plumed_assert(coeffs_pntr_in->numberOfDimensions()==nargs);
  plumed_assert(basisf_pntrs_in.size()==nargs);

  std::vector<double> args_values_trsfrm(nargs);
  std::vector< std::vector <double> > bf_values;
  //
  for(unsigned int k=0;k<nargs;k++){
    std::vector<double> tmp_val(basisf_pntrs_in[k]->getNumberOfBasisFunctions());
    std::vector<double> tmp_der(tmp_val.size());
    bool inside=true;
    basisf_pntrs_in[k]->getAllValues(args_values[k],args_values_trsfrm[k],inside,tmp_val,tmp_der);
    bf_values.push_back(tmp_val);
  }
  //
  size_t stride=1;
  size_t rank=0;
  if(comm_in!=NULL)
  {
    stride=comm_in->Get_size();
    rank=comm_in->Get_rank();
  }
  // loop over basis set
  for(size_t i=rank;i<coeffs_pntr_in->numberOfCoeffs();i+=stride){
    std::vector<unsigned int> indices=coeffs_pntr_in->getIndices(i);
    double bf_curr=1.0;
    for(unsigned int k=0;k<nargs;k++){
      bf_curr*=bf_values[k][indices[k]];
    }
    basisset_values[i] = bf_curr;
  }
  //
  if(comm_in!=NULL){
    comm_in->Sum(basisset_values);
  }
}


void LinearBasisSetExpansion::getBasisSetValues(const std::vector<double>& args_values, std::vector<double>& basisset_values) {
  getBasisSetValues(args_values,basisset_values,basisf_pntrs_, bias_coeffs_pntr_, &mycomm_);
}


void LinearBasisSetExpansion::setupUniformTargetDistribution() {
  std::vector<TargetDistribution*> targetdist_dummy(nargs_,NULL);
  setupSeperableTargetDistribution(targetdist_dummy);
}


void LinearBasisSetExpansion::setupTargetDistribution(const std::vector<std::string>& targetdist_keywords) {
  if(targetdist_keywords.size()!=1 && targetdist_keywords.size()!=nargs_){
    plumed_merror("the number of target distribution keywords needs to be either 1 or equal to the number of arguments");
  }
  std::vector<TargetDistribution*> targetdist_pntrs(targetdist_keywords.size(),NULL);
  for(unsigned int k=0;k<targetdist_keywords.size();k++){
    std::vector<std::string> words = Tools::getWords(targetdist_keywords[k]);
    if(words[0]=="UNIFORM"){
      targetdist_pntrs[k] = NULL;
    }
    else{
      targetdist_pntrs[k] = targetDistributionRegister().create(words);
    }
  }
  setupTargetDistribution(targetdist_pntrs);
  for(unsigned int k=0;k<targetdist_pntrs.size();k++){
    if(targetdist_pntrs[k]!=NULL){
      delete targetdist_pntrs[k];
    }
  }
}


void LinearBasisSetExpansion::setupTargetDistribution(const std::vector<TargetDistribution*>& targetdist_pntrs) {
  if(targetdist_pntrs.size()!=1 && targetdist_pntrs.size()!=nargs_){
    plumed_merror("the number of target distribution pointers needs to be either 1 or equal to the number of arguments");
  }
  if(targetdist_pntrs.size()==1){
    setupNonSeperableTargetDistribution(targetdist_pntrs[0]);
  }
  else{
    setupSeperableTargetDistribution(targetdist_pntrs);
  }
}


void LinearBasisSetExpansion::setupSeperableTargetDistribution(const std::vector<TargetDistribution*>& targetdist_pntrs) {
  plumed_massert(targetdist_pntrs.size()==nargs_,"number of target distribution does not match the number of basis functions");
  //
  std::vector< std::vector <double> > bf_integrals;
  for(unsigned int k=0;k<nargs_;k++){
    if(basisf_pntrs_[k]!=NULL){
      bf_integrals.push_back(basisf_pntrs_[k]->getTargetDistributionIntegrals(targetdist_pntrs[k]));
    }
    else{
      bf_integrals.push_back(basisf_pntrs_[k]->getUniformIntegrals());
    }
  }
  //
  for(size_t i=0;i<ncoeffs_;i++){
    std::vector<unsigned int> indices=bias_coeffs_pntr_->getIndices(i);
    double value = 1.0;
    for(unsigned int k=0;k<nargs_;k++){
      value*=bf_integrals[k][indices[k]];
    }
    coeffderivs_aver_ps_pntr_->setValue(i,value);
  }
}


void LinearBasisSetExpansion::setupNonSeperableTargetDistribution(const TargetDistribution* targetdist_pntr) {
  if(targetdist_pntr==NULL){
    setupUniformTargetDistribution();
    return;
  }
}


void LinearBasisSetExpansion::setupWellTemperedTargetDistribution(const double biasf, const std::vector<unsigned int>& nbins) {
  plumed_massert(biasf>1.0,"the value of the bias factor doesn't make sense, it should be larger than 1.0");
  biasf_=biasf;
  invbiasf_ = 1.0/biasf_;
  fes_wt_coeffs_pntr_ = new CoeffsVector(*bias_coeffs_pntr_);
  std::string fes_wt_label = bias_coeffs_pntr_->getLabel();
  if(fes_wt_label.find("coeffs")!=std::string::npos){
    fes_wt_label.replace(fes_wt_label.find("coeffs"), std::string("coeffs").length(), "fes_wt_coeffs");
  }
  else {
    fes_wt_label += "_fes_wt";
  }
  fes_wt_coeffs_pntr_->setLabels(fes_wt_label);
}


void LinearBasisSetExpansion::updateWellTemperedFESCoeffs() {
  plumed_assert(biasf_>1.0);
  FesWTCoeffs() = -BiasCoeffs() + invbiasf_*FesWTCoeffs();
}


}
}
