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
#include "tools/File.h"


#include <iostream>


namespace PLMD{
namespace bias{

VesBias::VesBias(const ActionOptions&ao):
Action(ao),
Bias(ao),
ncoeffssets_(0),
coeffs_pntrs_(0),
coeffderivs_aver_ps_pntrs_(0),
gradient_pntrs_(0),
hessian_pntrs_(0),
coeffderivs_aver_sampled(0),
coeffderivs_cov_sampled(0),
use_multiple_coeffssets_(false),
coeffs_fnames(0),
ncoeffs_total_(0),
optimizer_pntr_(NULL),
optimize_coeffs_(false),
compute_hessian_(false),
diagonal_hessian_(true),
targetdist_keywords_(0),
targetdist_pntrs_(0),
dynamic_targetdist_(false),
fname_coeffderivs_aver_ps("targetdist_averages.data"),
aver_counter(0.0),
kbt_(0.0),
welltemp_biasf_(1.0),
welltemp_targetdist_(false),
grid_bins_(0),
grid_min_(0),
grid_max_(0)
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

  if(keywords.exists("COEFFS")){
    parseVector("COEFFS",coeffs_fnames);
  }

  if(keywords.exists("TARGETDISTRIBUTION_AVERAGES_FILE")){
    parse("TARGETDISTRIBUTION_AVERAGES_FILE",fname_coeffderivs_aver_ps);
  }

  if(keywords.exists("GRID_BINS")){
    parseMultipleValues<unsigned int>("GRID_BINS",grid_bins_,getNumberOfArguments(),100);
  }

  if(keywords.exists("GRID_MIN") && keywords.exists("GRID_MAX")){
    parseMultipleValues("GRID_MIN",grid_min_,getNumberOfArguments());
    parseMultipleValues("GRID_MAX",grid_max_,getNumberOfArguments());
  }

  if(keywords.exists("TARGET_DISTRIBUTION")){
    std::string str_ps="";
    for(int i=1;;i++){
      if(!parseNumbered("TARGET_DISTRIBUTION",i,str_ps)){break;}
      targetdist_keywords_.push_back(str_ps);
    }
    str_ps="";
    parse("TARGET_DISTRIBUTION",str_ps);
    if(str_ps.size()>0){
      if(targetdist_keywords_.size()>0){
        plumed_merror("Either give a single target distribution using the TARGET_DISTRIBUTION keyword or multiple using numbered TARGET_DISTRIBUTION1,  TARGET_DISTRIBUTION2 keywords");
      }
      targetdist_keywords_.push_back(str_ps);
    }
  }

  if(keywords.exists("BIAS_FACTOR")){
    parse("BIAS_FACTOR",welltemp_biasf_);
    if(welltemp_biasf_<1.0){
      plumed_merror("the well-tempered bias factor doesn't make sense, it should be larger than 1.0");
    }
    if(welltemp_biasf_>1.0){
      welltemp_targetdist_=true;
      enableDynamicTargetDistribution();
      if(targetdist_keywords_.size()>0){
        plumed_merror("you cannot both specify a bias factor for a well-tempered target distribution using the BIAS_FACTOR keyword and give a target distribution with the TARGET_DISTRIBUTION keyword.");
      }
    }
  }


}


VesBias::~VesBias(){
  for(unsigned int i=0; i<coeffs_pntrs_.size(); i++){
    delete coeffs_pntrs_[i];
  }
  for(unsigned int i=0; i<coeffderivs_aver_ps_pntrs_.size(); i++){
    delete coeffderivs_aver_ps_pntrs_[i];
  }
  for(unsigned int i=0; i<gradient_pntrs_.size(); i++){
    delete gradient_pntrs_[i];
  }
  for(unsigned int i=0; i<hessian_pntrs_.size(); i++){
    delete hessian_pntrs_[i];
  }
}


void VesBias::registerKeywords( Keywords& keys ) {
  Bias::registerKeywords(keys);
  keys.add("optional","TEMP","the system temperature - this is needed if the MD code does not pass the temperature to PLUMED.");
  //
  keys.addOutputComponent("bias","default","the instantaneous value of the bias potential.");
  keys.addOutputComponent("force2","default","the instantaneous value of the squared force due to this bias potential.");
  //
  keys.reserve("optional","COEFFS","read-in the coefficents from files.");
  //
  keys.add("optional","TARGETDISTRIBUTION_AVERAGES_FILE","file for writing out the averages over the target distribution.");
  //
  //
  keys.reserve("numbered","TARGET_DISTRIBUTION","the target distribution to be used.");
    //
  keys.reserve("optional","BIAS_FACTOR","the bias factor to be used for the well-tempered target distribution.");
  //
  keys.reserve("optional","GRID_BINS","the number of bins used for the grid. The default value is 100 bins per dimension.");
  keys.reserve("optional","GRID_MIN","the lower bounds used for the grid.");
  keys.reserve("optional","GRID_MAX","the upper bounds used for the grid.");
}


void VesBias::useInitialCoeffsKeywords(Keywords& keys) {
  keys.use("COEFFS");
}


void VesBias::useTargetDistributionKeywords(Keywords& keys) {
  keys.use("TARGET_DISTRIBUTION");
}


void VesBias::useWellTemperdKeywords(Keywords& keys) {
  keys.use("BIAS_FACTOR");
}


void VesBias::useGridBinKeywords(Keywords& keys) {
  keys.use("GRID_BINS");
}


void VesBias::useGridLimitsKeywords(Keywords& keys) {
  keys.use("GRID_MIN");
  keys.use("GRID_MAX");
}


void VesBias::addCoeffsSet(const std::vector<std::string>& dimension_labels,const std::vector<unsigned int>& indices_shape) {
  CoeffsVector* coeffs_pntr_tmp = new CoeffsVector("coeffs",dimension_labels,indices_shape,comm,true);
  initializeCoeffs(coeffs_pntr_tmp);
}


void VesBias::addCoeffsSet(std::vector<Value*>& args,std::vector<BasisFunctions*>& basisf) {
  CoeffsVector* coeffs_pntr_tmp = new CoeffsVector("coeffs",args,basisf,comm,true);
  initializeCoeffs(coeffs_pntr_tmp);
}


void VesBias::addCoeffsSet(CoeffsVector* coeffs_pntr_in) {
  initializeCoeffs(coeffs_pntr_in);
}


void VesBias::initializeCoeffs(CoeffsVector* coeffs_pntr_in) {
  //
  coeffs_pntr_in->linkVesBias(this);
  //
  std::string label;
  if(!use_multiple_coeffssets_ && ncoeffssets_==1){
    plumed_merror("you are not allowed to use multiple coefficient sets");
  }
  //
  label = getCoeffsSetLabelString("coeffs",ncoeffssets_);
  coeffs_pntr_in->setLabels(label);

  coeffs_pntrs_.push_back(coeffs_pntr_in);
  CoeffsVector* aver_ps_tmp = new CoeffsVector(*coeffs_pntr_in);
  label = getCoeffsSetLabelString("ps-aver",ncoeffssets_);
  aver_ps_tmp->setLabels(label);
  aver_ps_tmp->setValues(0.0);
  coeffderivs_aver_ps_pntrs_.push_back(aver_ps_tmp);
  //
  CoeffsVector* gradient_tmp = new CoeffsVector(*coeffs_pntr_in);
  label = getCoeffsSetLabelString("gradient",ncoeffssets_);
  gradient_tmp->setLabels(label);
  gradient_pntrs_.push_back(gradient_tmp);
  //
  label = getCoeffsSetLabelString("hessian",ncoeffssets_);
  CoeffsMatrix* hessian_tmp = new CoeffsMatrix(label,coeffs_pntr_in,comm,diagonal_hessian_);
  hessian_pntrs_.push_back(hessian_tmp);
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


void VesBias::clearCoeffsPntrsVector() {
  coeffs_pntrs_.clear();
}


void VesBias::readCoeffsFromFiles() {
  plumed_assert(ncoeffssets_>0);
  plumed_massert(keywords.exists("COEFFS"),"you are not allowed to use this function as the COEFFS keyword is not enabled");
  if(coeffs_fnames.size()>0){
    plumed_massert(coeffs_fnames.size()==ncoeffssets_,"COEFFS keyword is of the wrong size");
    if(ncoeffssets_==1){
      log.printf("  Read in coefficents from file ");
    }
    else{
      log.printf("  Read in coefficents from files:\n");
    }
    for(unsigned int i=0; i<ncoeffssets_; i++){
      IFile ifile;
      ifile.link(*this);
      ifile.open(coeffs_fnames[i]);
      if(!ifile.FieldExist(coeffs_pntrs_[i]->getDataLabel())){
        std::string error_msg = "Problem with reading coefficents from file " + ifile.getPath() + ": no field with name " + coeffs_pntrs_[i]->getDataLabel() + "\n";
        plumed_merror(error_msg);
      }
      size_t ncoeffs_read = coeffs_pntrs_[i]->readFromFile(ifile,false,false);
      coeffs_pntrs_[i]->setIterationCounterAndTime(0,getTime());
      if(ncoeffssets_==1){
        log.printf("%s (read %zu of %zu values)\n", ifile.getPath().c_str(),ncoeffs_read,coeffs_pntrs_[i]->numberOfCoeffs());
      }
      else{
        log.printf("   coefficent %u: %s (read %zu of %zu values)\n",i,ifile.getPath().c_str(),ncoeffs_read,coeffs_pntrs_[i]->numberOfCoeffs());
      }
      ifile.close();
    }
  }
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


void VesBias::setCoeffsDerivsOverTargetDist(const CoeffsVector& coeffderivs_aver_ps, const unsigned coeffs_id) {
  CoeffDerivsAverTargetDist(coeffs_id) = coeffderivs_aver_ps;
}


void VesBias::setCoeffsDerivsOverTargetDistToZero(const unsigned coeffs_id) {
  CoeffDerivsAverTargetDist(coeffs_id).setAllValuesToZero();
}


void VesBias::checkThatTemperatureIsGiven() {
  if(kbt_==0.0){
    std::string err_msg = "VES bias " + getLabel() + " of type " + getName() + ": the temperature is needed so you need to give it using the TEMP keyword as the MD engine does not pass it to PLUMED.";
    plumed_merror(err_msg);
  }
}


void VesBias::linkOptimizer(Optimizer* optimizer_pntr_in) {
  //
  if(optimizer_pntr_==NULL){
    optimizer_pntr_ = optimizer_pntr_in;
  }
  else {
    std::string err_msg = "VES bias " + getLabel() + " of type " + getName() + " has already been linked with optimizer " + optimizer_pntr_->getLabel() + " of type " + optimizer_pntr_->getName() + ". You cannot link two optimizer to the same VES bias.";
    plumed_merror(err_msg);
  }
  checkThatTemperatureIsGiven();
  optimize_coeffs_ = true;
}


void VesBias::enableHessian(const bool diagonal_hessian) {
  compute_hessian_=true;
  diagonal_hessian_=diagonal_hessian;
  coeffderivs_cov_sampled.clear();
  for (unsigned int i=0; i<ncoeffssets_; i++){
    delete hessian_pntrs_[i];
    std::string label = getCoeffsSetLabelString("hessian",i);
    hessian_pntrs_[i] = new CoeffsMatrix(label,coeffs_pntrs_[i],comm,diagonal_hessian_);
    std::vector<double> cov_sampled_tmp;
    cov_sampled_tmp.assign(hessian_pntrs_[i]->getSize(),0.0);
    coeffderivs_cov_sampled.push_back(cov_sampled_tmp);
  }
}


void VesBias::disableHessian() {
  compute_hessian_=false;
  diagonal_hessian_=true;
  coeffderivs_cov_sampled.clear();
  for (unsigned int i=0; i<ncoeffssets_; i++){
    delete hessian_pntrs_[i];
    std::string label = getCoeffsSetLabelString("hessian",i);
    hessian_pntrs_[i] = new CoeffsMatrix(label,coeffs_pntrs_[i],comm,diagonal_hessian_);
    std::vector<double> cov_sampled_tmp;
    cov_sampled_tmp.assign(hessian_pntrs_[i]->getSize(),0.0);
    coeffderivs_cov_sampled.push_back(cov_sampled_tmp);
  }
}


void VesBias::apply() {
  Bias::apply();
  aver_counter += 1.0;
}


std::string VesBias::getCoeffsSetLabelString(const std::string& type, const unsigned int coeffs_id) {
  std::string label_prefix = getLabel() + ".";
  std::string label_postfix = "";
  if(use_multiple_coeffssets_){
    Tools::convert(coeffs_id,label_postfix);
    label_postfix = "-" + label_postfix;
  }
  return label_prefix+type+label_postfix;
}


void VesBias::updateTargetDistributions() {}


void VesBias::writeCoeffDerivsAverTargetDistToFile(const bool append, const unsigned int iteration) {
  getCoeffDerivsAverTargetDistPntr()->setIterationCounterAndTime(iteration,this->getTime());
  getCoeffDerivsAverTargetDistPntr()->writeToFile(fname_coeffderivs_aver_ps,true,append,static_cast<Action*>(this));
}


void VesBias::setGridBins(const std::vector<unsigned int>& grid_bins_in) {
  plumed_massert(grid_bins_in.size()==getNumberOfArguments(),"the number of grid bins given doesn't match the number of arguments");
  grid_bins_=grid_bins_in;
}


void VesBias::setGridBins(const unsigned int nbins) {
  std::vector<unsigned int> grid_bins_in(getNumberOfArguments(),nbins);
  grid_bins_=grid_bins_in;
}


void VesBias::setGridMin(const std::vector<double>& grid_min_in) {
  plumed_massert(grid_min_in.size()==getNumberOfArguments(),"the number of lower bounds given for the grid doesn't match the number of arguments");
  grid_min_=grid_min_in;
}


void VesBias::setGridMax(const std::vector<double>& grid_max_in) {
  plumed_massert(grid_max_in.size()==getNumberOfArguments(),"the number of upper bounds given for the grid doesn't match the number of arguments");
  grid_max_=grid_max_in;
}


}
}
