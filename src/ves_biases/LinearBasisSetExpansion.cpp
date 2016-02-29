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
  const double beta_in,
  Communicator& cc,
  std::vector<Value*> args_pntrs_in,
  std::vector<BasisFunctions*> basisf_pntrs_in,
  CoeffsVector* bias_coeffs_pntr_in):
label_(label),
action_pntr_(NULL),
vesbias_pntr_(NULL),
mycomm_(cc),
serial_(false),
beta_(beta_in),
args_pntrs_(args_pntrs_in),
nargs_(args_pntrs_.size()),
basisf_pntrs_(basisf_pntrs_in),
nbasisf_(basisf_pntrs_.size()),
bias_coeffs_pntr_(bias_coeffs_pntr_in),
ncoeffs_(0),
coeffderivs_aver_ps_pntr_(NULL),
fes_wt_coeffs_pntr_(NULL),
welltemp_biasf_(-1.0),
inv_welltemp_biasf_(-1.0),
beta_prime_(0.0),
grid_bins_(nargs_,100),
bias_grid_pntr_(NULL),
fes_grid_pntr_(NULL),
ps_grid_pntr_(NULL),
welltemp_ps_grid_pntr_(NULL)
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
  if(coeffderivs_aver_ps_pntr_!=NULL){
    delete coeffderivs_aver_ps_pntr_;
  }
  if(fes_wt_coeffs_pntr_!=NULL){
    delete fes_wt_coeffs_pntr_;
  }
  if(bias_grid_pntr_!=NULL){
    delete bias_grid_pntr_;
  }
  if(fes_grid_pntr_!=NULL){
    delete fes_grid_pntr_;
  }
  if(ps_grid_pntr_!=NULL){
    delete ps_grid_pntr_;
  }
  if(welltemp_ps_grid_pntr_!=NULL){
    delete welltemp_ps_grid_pntr_;
  }

}


void LinearBasisSetExpansion::linkVesBias(bias::VesBias* vesbias_pntr_in) {
  vesbias_pntr_ = vesbias_pntr_in;
  action_pntr_ = static_cast<Action*>(vesbias_pntr_in);
}


void LinearBasisSetExpansion::linkAction(Action* action_pntr_in) {
  action_pntr_ = action_pntr_in;
}


void LinearBasisSetExpansion::setGridBins(const std::vector<unsigned int>& grid_bins_in) {
  plumed_massert(grid_bins_in.size()==nargs_,"the number of grid bins given doesn't match the number of arguments");
  grid_bins_=grid_bins_in;
}


void LinearBasisSetExpansion::setGridBins(const unsigned int nbins) {
  std::vector<unsigned int> grid_bins_in(nargs_,nbins);
  grid_bins_=grid_bins_in;
}


Grid* LinearBasisSetExpansion::setupGeneralGrid(const std::string label_suffix, const std::vector<unsigned int>& nbins, const bool usederiv) {
  plumed_assert(nbins.size()==nargs_);
  std::vector<std::string> min(nargs_);
  std::vector<std::string> max(nargs_);
  for(unsigned int k=0;k<nargs_;k++){
    Tools::convert(basisf_pntrs_[k]->intervalMin(),min[k]);
    Tools::convert(basisf_pntrs_[k]->intervalMax(),max[k]);
  }
  Grid* grid_pntr = new Grid(label_+"."+label_suffix,args_pntrs_,min,max,nbins,false,usederiv);
  return grid_pntr;
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
  plumed_assert(coeffsderivs_values.size()==coeffs_pntr_in->numberOfCoeffs());

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
  if(nargs_>1 && targetdist_pntrs.size()==1){
    setupNonSeperableTargetDistribution(targetdist_pntrs[0]);
  }
  else{
    setupSeperableTargetDistribution(targetdist_pntrs);
  }
}


void LinearBasisSetExpansion::setupSeperableTargetDistribution(const std::vector<TargetDistribution*>& targetdist_pntrs) {
  plumed_massert(targetdist_pntrs.size()==nargs_,"number of target distribution does not match the number of basis functions");
  //
  for(unsigned int k=0;k<nargs_;k++){
    std::vector<std::string> min(1);
    std::vector<std::string> max(1);
    Tools::convert(basisf_pntrs_[k]->intervalMin(),min[0]);
    Tools::convert(basisf_pntrs_[k]->intervalMax(),max[0]);
    std::vector<unsigned int> nbins(1);
    nbins[0]=300;
    std::string ks; Tools::convert(k+1,ks);
    std::string filename = "targetdist-" + ks + ".data";
    if(targetdist_pntrs[k]!=NULL){
      targetdist_pntrs[k]->writeDistributionToFile(filename,min,max,nbins);
    }

  }
  //
  std::vector< std::vector <double> > bf_integrals(0);
  for(unsigned int k=0;k<nargs_;k++){
    if(targetdist_pntrs[k]!=NULL){
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
  Grid* ps_grid_pntr = setupGeneralGrid("ps",grid_bins_,false);
  targetdist_pntr->calculateDistributionOnGrid(ps_grid_pntr);
  calculateCoeffDerivsAverFromGrid(ps_grid_pntr);
  TargetDistribution::writeProbGridToFile("targetdist.data",ps_grid_pntr,true);
  delete ps_grid_pntr;
}


void LinearBasisSetExpansion::setupWellTemperedTargetDistribution(const double biasf) {
  plumed_massert(biasf>1.0,"the value of the bias factor doesn't make sense, it should be larger than 1.0");
  welltemp_biasf_=biasf;
  inv_welltemp_biasf_ = 1.0/welltemp_biasf_;
  beta_prime_ = beta_/welltemp_biasf_;
  fes_wt_coeffs_pntr_ = new CoeffsVector(*bias_coeffs_pntr_);
  std::string fes_wt_label = bias_coeffs_pntr_->getLabel();
  if(fes_wt_label.find("coeffs")!=std::string::npos){
    fes_wt_label.replace(fes_wt_label.find("coeffs"), std::string("coeffs").length(), "fes_wt_coeffs");
  }
  else {
    fes_wt_label += "_fes_wt";
  }
  fes_wt_coeffs_pntr_->setLabels(fes_wt_label);
  welltemp_ps_grid_pntr_ = setupGeneralGrid("ps_wt",grid_bins_,false);
}


void LinearBasisSetExpansion::updateWellTemperedPsGrid() {
  double norm = 0.0;
  size_t stride=mycomm_.Get_size();
  size_t rank=mycomm_.Get_rank();
  for(unsigned int l=rank; l<welltemp_ps_grid_pntr_->getSize(); l+=stride){
    std::vector<double> args_values = welltemp_ps_grid_pntr_->getPoint(l);
    double value = -beta_prime_ * getFES_WellTempered(args_values,false);
    value = exp(value);
    norm += value;
    welltemp_ps_grid_pntr_->setValue(l,value);
  }
  norm = 1.0/(welltemp_ps_grid_pntr_->getBinVolume()*norm);
  welltemp_ps_grid_pntr_->scaleAllValuesAndDerivatives(norm);
  welltemp_ps_grid_pntr_->mpiSumValuesAndDerivatives(mycomm_);
}


void LinearBasisSetExpansion::updateWellTemperedTargetDistribution() {
  plumed_assert(welltemp_biasf_>1.0);
  // Update FES WT coeffs
  FesWTCoeffs() = -BiasCoeffs() + inv_welltemp_biasf_*FesWTCoeffs();
  // Update p(s) grid
  updateWellTemperedPsGrid();
  // calcuate coeffs derivs from grid
  calculateCoeffDerivsAverFromGrid(welltemp_ps_grid_pntr_);
}


void LinearBasisSetExpansion::calculateCoeffDerivsAverFromGrid(const Grid* ps_grid_pntr, const bool normalize_dist) {
  plumed_assert(ps_grid_pntr!=NULL);
  std::vector<double> coeffderivs_aver_ps(ncoeffs_,0.0);
  double sum_grid = 0.0;
  double binVol = ps_grid_pntr->getBinVolume();
  size_t stride=mycomm_.Get_size();
  size_t rank=mycomm_.Get_rank();
  for(unsigned int l=rank; l<ps_grid_pntr->getSize(); l+=stride){
    std::vector<double> args_values = ps_grid_pntr->getPoint(l);
    std::vector<double> basisset_values(ncoeffs_);
    getBasisSetValues(args_values,basisset_values,false);
    double weight = ps_grid_pntr->getValue(l)*binVol;
    sum_grid += weight;
    for(unsigned int i=0; i<ncoeffs_; i++){
      coeffderivs_aver_ps[i] += weight*basisset_values[i];
    }
  }
  if(normalize_dist){
    mycomm_.Sum(sum_grid);
    for(unsigned int i=rank; i<ncoeffs_; i+=stride){
      coeffderivs_aver_ps[i] /= sum_grid;
    }
  }
  // the overall constant;
  mycomm_.Sum(coeffderivs_aver_ps);
  coeffderivs_aver_ps[0] = 1.0;
  CoeffDerivsAverTargetDist() = coeffderivs_aver_ps;
}

}

}
