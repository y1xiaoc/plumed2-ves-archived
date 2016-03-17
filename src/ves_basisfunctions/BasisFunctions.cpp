/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2011-2014 The plumed team
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
#include "BasisFunctions.h"
#include "ves_targetdistributions/TargetDistribution.h"
#include "ves_targetdistributions/TargetDistributionRegister.h"
#include "ves_biases/VesBias.h"


namespace PLMD{

void BasisFunctions::registerKeywords(Keywords& keys){
  Action::registerKeywords(keys);
  ActionWithValue::registerKeywords(keys);
  keys.add("compulsory","ORDER","The order of the basis functions.");
  keys.add("compulsory","INTERVAL_MIN","the minimum of the interval on which the basis functions are defined");
  keys.add("compulsory","INTERVAL_MAX","the maximum of the interval on which the basis functions are defined");
  keys.add("optional","NGRID_POINTS","the number of grid points used for numerical integrals");
  keys.addFlag("DEBUG_INFO",false,"print out more detailed information about the basis set, useful for debugging");
  keys.addFlag("NUMERICAL_INTEGRALS",false,"calculate basis function integral for the uniform distribution numerically");
}


BasisFunctions::BasisFunctions(const ActionOptions&ao):
Action(ao),
ActionWithValue(ao),
print_debug_info_(false),
has_been_set(false),
description_("Undefined"),
type_("Undefined"),
norder_(0),
nbasis_(1),
bf_label_prefix_("f"),
bf_labels_(nbasis_,"f0"),
periodic_(false),
interval_bounded_(true),
interval_intrinsic_min_(1.0),
interval_intrinsic_max_(-1.0),
interval_intrinsic_range_(0.0),
interval_intrinsic_mean_(0.0),
interval_min_(0.0),
interval_max_(0.0),
interval_range_(0.0),
interval_mean_(0.0),
argT_derivf_(1.0),
numerical_uniform_integrals_(false),
nbins_(1001),
uniform_integrals_(nbasis_,0.0),
vesbias_pntr_(NULL),
action_pntr_(NULL)
{
  bf_keywords_.push_back(getName());
  parse("ORDER",norder_); addKeywordToList("ORDER",norder_);
  nbasis_=norder_+1;
  //
  std::string str_imin_; std::string str_imax_;
  parse("INTERVAL_MIN",str_imin_); addKeywordToList("INTERVAL_MIN",str_imin_);
  parse("INTERVAL_MAX",str_imax_); addKeywordToList("INTERVAL_MAX",str_imax_);
  Tools::convert(str_imin_,interval_min_); Tools::convert(str_imax_,interval_max_);
  if(interval_min_>interval_max_){plumed_merror("INTERVAL_MIN and INTERVAL_MIX are not correctly defined");}
  //
  parseFlag("DEBUG_INFO",print_debug_info_);
  parseFlag("NUMERICAL_INTEGRALS",numerical_uniform_integrals_);
  parse("NGRID_POINTS",nbins_);
  // log.printf(" %s \n",getKeywordString().c_str());

}


void BasisFunctions::setNumberOfBasisFunctions(const unsigned int nbasis_in) {
  nbasis_=nbasis_in;
  bf_labels_.assign(nbasis_,"");
  uniform_integrals_.assign(nbasis_,0.0);
}


void BasisFunctions::setupInterval() {
  // if(!intervalBounded()){plumed_merror("setupInterval() only works for bounded interval");}
  interval_intrinsic_range_ = interval_intrinsic_max_-interval_intrinsic_min_;
  interval_intrinsic_mean_  = 0.5*(interval_intrinsic_max_+interval_intrinsic_min_);
  interval_range_ = interval_max_-interval_min_;
  interval_mean_  = 0.5*(interval_max_+interval_min_);
  argT_derivf_ = interval_intrinsic_range_/interval_range_;
}


void BasisFunctions::setupLabels() {
  for(unsigned int i=0; i < nbasis_;i++){
    std::string is; Tools::convert(i,is);
    bf_labels_[i]=bf_label_prefix_+is+"(s)";
  }
}


void BasisFunctions::setupUniformIntegrals() {
  numerical_uniform_integrals_=true;
  numericalUniformIntegrals();
}


void BasisFunctions::setupBF() {
  checkRead();
  if(interval_intrinsic_min_>interval_intrinsic_max_){plumed_merror("setupBF: default intervals are not correctly set");}
  setupInterval();
  setupLabels();
  if(bf_labels_.size()==1){plumed_merror("setupBF: the labels of the basis functions are not correct.");}
  if(!numerical_uniform_integrals_){setupUniformIntegrals();}
  else{numericalUniformIntegrals();}
  if(uniform_integrals_.size()==1){plumed_merror("setupBF: the integrals of the basis functions is not correct.");}
  if(type_=="Undefined"){plumed_merror("setupBF: the type of the basis function is not defined.");}
  if(description_=="Undefined"){plumed_merror("setupBF: the description of the basis function is not defined.");}
  has_been_set=true;
  printInfo();
}


void BasisFunctions::printInfo() const {
  if(!has_been_set){plumed_merror("the basis set has not be setup correctly");}
  log.printf("  One-dimensional basis set\n");
  log.printf("   Description: %s\n",description_.c_str());
  log.printf("   Type: %s\n",type_.c_str());
  if(periodic_){log.printf("   The basis functions are periodic\n");}
  log.printf("   Order of basis set: %u\n",norder_);
  log.printf("   Number of basis functions: %u\n",nbasis_);
  log.printf("   Interval of basis set: %f to %f\n",interval_min_,interval_max_);
  log.printf("   Description of basis functions:\n");
  for(unsigned int i=0; i < nbasis_;i++){log.printf("    %2u       %10s\n",i,bf_labels_[i].c_str());}
  //
  if(print_debug_info_){
    log.printf("  Debug information:\n");
    log.printf("   Default interval of basis set: [%f,%f]\n",interval_intrinsic_min_,interval_intrinsic_max_);
    log.printf("   Default interval of basis set: range=%f,  mean=%f\n",interval_intrinsic_range_,interval_intrinsic_mean_);
    log.printf("   Defined interval of basis set: [%f,%f]\n",interval_min_,interval_max_);
    log.printf("   Defined interval of basis set: range=%f,  mean=%f\n",interval_range_,interval_mean_);
    log.printf("   Derivative factor due to interval translation: %f\n",argT_derivf_);
    log.printf("   Integral of basis functions over the interval:\n");
    if(numerical_uniform_integrals_){log.printf("   Note: calculated numerically\n");}
    for(unsigned int i=0; i < nbasis_;i++){log.printf("    %2u       %16.10f\n",i,uniform_integrals_[i]);}
    log.printf("   --------------------------\n");
  }
}


void BasisFunctions::linkVesBias(bias::VesBias* vesbias_pntr_in){
  vesbias_pntr_ = vesbias_pntr_in;
  action_pntr_ = static_cast<Action*>(vesbias_pntr_in);
}


void BasisFunctions::linkAction(Action* action_pntr_in){
  action_pntr_ = action_pntr_in;
}


void BasisFunctions::numericalUniformIntegrals() {
  double h=(interval_max_-interval_min_)/nbins_;
  //
  bool dummy_bool=true;
  double dummy_dbl=0.0;
  for(unsigned int i=0; i < nbasis_;i++){
    // Trapezoidal rule on a uniform grid with Nbins+1 grid points
    double sum=0.0;
    for(unsigned int k=0; k < nbins_;k++){
      double x1 = interval_min_+(k)*h;
      double x2 = interval_min_+(k+1)*h;
      double v1 = getValue(x1,i,dummy_dbl,dummy_bool);
      double v2 = getValue(x2,i,dummy_dbl,dummy_bool);
      sum = sum + (v1+v2);
    }
    // norm with the "volume of the interval"
    // assume that the first function is the constant
    uniform_integrals_[0] = getValue(0.0,0,dummy_dbl,dummy_bool);
    uniform_integrals_[i] = (0.5*h*sum)/interval_range_;
  }
}


std::vector<double> BasisFunctions::numericalTargetDistributionIntegrals(const TargetDistribution* targetdist_pntr) const {
  if(targetdist_pntr==NULL){
    return getUniformIntegrals();
  }
  plumed_massert(targetdist_pntr->getDimension()==1,"the target distribution should be one dimensional");
  //
  double h=(interval_max_-interval_min_)/nbins_;
  std::vector<double> targetdist_integrals(nbasis_,0.0);
  //
  std::vector<double> weights(nbins_+1,0.0);
  for(unsigned int k=0; k < (nbins_+1); k++){
    std::vector<double> x1(1);
    x1[0] = interval_min_+(k)*h;
    weights[k] = targetdist_pntr->getValue(x1);
  }
  if( !(targetdist_pntr->isNormalized()) ){
    double norm = 0.0;
    for(unsigned int k=0; k < nbins_; k++){
      norm = norm + (weights[k] + weights[k+1]);
    }
    norm = 1.0/(0.5*h*norm);
    for(unsigned int k=0; k < (nbins_+1); k++){
      weights[k] = norm*weights[k];
    }
  }
  //
  bool dummy_bool=true;
  double dummy_dbl=0.0;
  for(unsigned int i=0; i < nbasis_;i++){
    // Trapezoidal rule on a uniform grid with Nbins+1 grid points
    double sum=0.0;
    for(unsigned int k=0; k < nbins_; k++){
      double x1 = interval_min_+(k)*h;
      double x2 = interval_min_+(k+1)*h;
      double v1 = weights[k]*getValue(x1,i,dummy_dbl,dummy_bool);
      double v2 = weights[k+1]*getValue(x2,i,dummy_dbl,dummy_bool);
      sum = sum + (v1+v2);
    }
    // norm with the "volume of the interval"
    targetdist_integrals[i] = (0.5*h*sum);
  }
  //
  // assume that the first function is the constant
  targetdist_integrals[0] = getValue(0.0,0,dummy_dbl,dummy_bool);
  return targetdist_integrals;

}


std::vector<double> BasisFunctions::getTargetDistributionIntegrals(const std::string& targetdist_keyword) const {
  std::vector<std::string> words = Tools::getWords(targetdist_keyword);
  TargetDistribution* targetdist_pntr = NULL;
  if(words[0]=="UNIFORM"){
    targetdist_pntr = NULL;
  }
  else{
    targetdist_pntr = targetDistributionRegister().create(words);
  }
  std::vector<double> integrals = getTargetDistributionIntegrals(targetdist_pntr);
  if(targetdist_pntr!=NULL){
    delete targetdist_pntr;
  }
  return integrals;
}


std::string BasisFunctions::getKeywordString() const {
  std::string str_keywords=bf_keywords_[0];
  for(unsigned int i=1; i<bf_keywords_.size();i++){str_keywords+=" "+bf_keywords_[i];}
  return str_keywords;
}


}
