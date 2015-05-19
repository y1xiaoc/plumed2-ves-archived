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


namespace PLMD{

BasisFunctions::BasisFunctions(const ActionOptions&ao):
Action(ao),
ActionWithValue(ao),
print_debug_info_(false),
has_been_set(false),
description_("Undefined"),
type_("Undefined"),
norder_(0),
nbasis_(1),
bf_description_prefix_("f"),
bf_description_(nbasis_,"f0"),
periodic_(false),
interval_bounded_(true),
interval_default_min_(1.0),
interval_default_max_(-1.0),
interval_default_range_(0.0),
interval_default_mean_(0.0),
interval_min_(0.0),
interval_max_(0.0),
interval_range_(0.0),
interval_mean_(0.0),
argT_derivf_(1.0),
numerical_bf_integrals_(false),
bf_integrals_(nbasis_,0.0)
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
 parseFlag("NUMERICAL_BF_INTEGRALS",numerical_bf_integrals_); 
 // log.printf(" %s \n",getKeywordString().c_str());
}

void BasisFunctions::registerKeywords( Keywords& keys ){
  Action::registerKeywords(keys);
  ActionWithValue::registerKeywords(keys);
  keys.add("compulsory","ORDER","The order of the basis functions.");
  keys.add("compulsory","INTERVAL_MIN","the minimum of the interval on which the basis functions are defined");
  keys.add("compulsory","INTERVAL_MAX","the maximum of the interval on which the basis functions are defined");
  keys.addFlag("DEBUG_INFO",false,"print out more detailed information about the basis set, useful for debugging");
  keys.addFlag("NUMERICAL_BF_INTEGRALS",false,"calculate basis function integral over the interval numerically");
}

void BasisFunctions::setupInterval(){
 // if(!intervalBounded()){plumed_merror("setupInterval() only works for bounded interval");}
 interval_default_range_ = interval_default_max_-interval_default_min_;
 interval_default_mean_  = 0.5*(interval_default_max_+interval_default_min_);
 interval_range_ = interval_max_-interval_min_;
 interval_mean_  = 0.5*(interval_max_+interval_min_);
 argT_derivf_ = interval_default_range_/interval_range_;
}

double BasisFunctions::translateArgument(const double arg, bool& inside_interval)
{
 inside_interval=true;
 double argT = (arg-interval_mean_)*argT_derivf_;
 if(argT < interval_default_min_){
  inside_interval=false;
  argT=interval_default_min_;
 }
 else if(argT > interval_default_max_){
  inside_interval=false;
  argT=interval_default_max_;
 }
 return argT;
}

void BasisFunctions::apply(){}

void BasisFunctions::calculate(){}

void BasisFunctions::setupDescription()
{
 bf_description_.resize(nbasis_);
 for(unsigned int i=0; i < nbasis_;i++)
 {
  std::string is; Tools::convert(i,is);
  bf_description_[i]=bf_description_prefix_+is+"(s)";
 }
}

void BasisFunctions::setupBFIntegrals()
{
 numerical_bf_integrals_=true;
 numericalBFIntegrals();
}

void BasisFunctions::setupBF()
{
 checkRead();
 if(interval_default_min_>interval_default_max_){plumed_merror("setupBF: default intervals are not correctly set");}
 setupInterval();
 setupDescription();
 if(bf_description_.size()==1){plumed_merror("setupBF: the description of the basis functions is not correct.");}
 if(!numerical_bf_integrals_){setupBFIntegrals();}
 else{numericalBFIntegrals();}
 if(bf_integrals_.size()==1){plumed_merror("setupBF: the integrals of the basis functions is not correct.");}
 if(type_=="Undefined"){plumed_merror("setupBF: the type of the basis function is not defined.");}
 if(description_=="Undefined"){plumed_merror("setupBF: the description of the basis function is not defined.");}
 has_been_set=true;
}

void BasisFunctions::printInfo()
{
 if(!has_been_set){plumed_merror("the basis set has not be setup correctly");}
 log.printf("  One-dimensional basis set\n");
 log.printf("   Description: %s\n",description_.c_str());
 log.printf("   Type: %s\n",type_.c_str());
 if(periodic_){log.printf("   The basis functions are periodic\n");}
 log.printf("   Order of basis set: %d\n",norder_);
 log.printf("   Number of basis functions: %d\n",nbasis_);
 log.printf("   Interval of basis set: %f to %f\n",interval_min_,interval_max_);
 log.printf("   Description of basis functions:\n");
 for(unsigned int i=0; i < nbasis_;i++){log.printf("    %2d       %10s\n",i,bf_description_[i].c_str());}
 
 if(print_debug_info_)
 {
  log.printf("  Debug information:\n");
  log.printf("   Default interval of basis set: [%f,%f]\n",interval_default_min_,interval_default_max_);
  log.printf("   Default interval of basis set: range=%f,  mean=%f\n",interval_default_range_,interval_default_mean_);
  log.printf("   Defined interval of basis set: [%f,%f]\n",interval_min_,interval_max_);
  log.printf("   Defined interval of basis set: range=%f,  mean=%f\n",interval_range_,interval_mean_);
  log.printf("   Derivative factor due to interval translation: %f\n",argT_derivf_);
  log.printf("   Integral of basis functions over the interval:\n");
  if(numerical_bf_integrals_){log.printf("   Note: calculated numerically\n");}
  for(unsigned int i=0; i < nbasis_;i++){log.printf("    %2d       %16.10f\n",i,bf_integrals_[i]);}
  log.printf("   --------------------------\n");
 }
}

void BasisFunctions::numericalBFIntegrals()
{
 unsigned int nbins=1001;
 double h=(interval_max_-interval_min_)/nbins;
 bf_integrals_.assign(nbasis_,0.0);

 bool dummy_bool=true;
 double dummy_dbl=0.0;
 for(unsigned int i=0; i < nbasis_;i++)
 {
  // Trapezoidal rule on a uniform grid with Nbins+1 grid points
  double sum=0.0;
  for(unsigned int k=0; k < nbins;k++)
  {
   double x1 = interval_min_+(k)*h;
   double x2 = interval_min_+(k+1)*h;
   double v1 = getValue(x1,i,dummy_dbl,dummy_bool);
   double v2 = getValue(x2,i,dummy_dbl,dummy_bool);
   sum = sum + (v1+v2);
  }
  // norm with the "volume of the interval"
  bf_integrals_[i] = (0.5*h*sum)/interval_range_; 
 } 
}

template<typename T>
void BasisFunctions::addKeywordToList(const std::string keyword, const T value)
{
 std::string str_value;
 Tools::convert(value,str_value);
 bf_keywords_.push_back(keyword+"="+str_value);
}

void BasisFunctions::addKeywordToList(const std::string keyword, const bool value)
{
 if(value){bf_keywords_.push_back(keyword);}
}

std::string BasisFunctions::getKeywordString()
{
 std::string str_keywords=bf_keywords_[0];
 for(unsigned int i=1; i<bf_keywords_.size();i++){str_keywords+=" "+bf_keywords_[i];}
 return str_keywords;
}

}


