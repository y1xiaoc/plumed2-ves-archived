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
namespace BasisFunctions{

BasisFunctions::BasisFunctions(const ActionOptions&ao):
Action(ao),
ActionWithValue(ao),
has_been_set(false),
type_("Undefined"),
norder_(0),
nbasis_(1),
bf_description_(nbasis_,"1"),
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
bf_integrals_(nbasis_,0.0)
{
 parse("ORDER",norder_);
 parse("INTERVAL_MIN",interval_min_);
 parse("INTERVAL_MAX",interval_max_);
 if(interval_min_>interval_max_){error("INTERVAL_MIN and INTERVAL_MIX are not correctly defined");}
}

void BasisFunctions::registerKeywords( Keywords& keys ){
  Action::registerKeywords(keys);
  ActionWithValue::registerKeywords(keys);
  keys.add("compulsory","ORDER","The order of the basis functions.");
  keys.add("compulsory","INTERVAL_MIN","the minimum of the interval on which the basis functions are defined");
  keys.add("compulsory","INTERVAL_MAX","the maximum of the interval on which the basis functions are defined");
}

void BasisFunctions::setupInterval(){
 // if(!intervalBounded()){error("setupInterval() only works for bounded interval");}
 interval_default_range_ = interval_default_max_-interval_default_min_;
 interval_default_mean_  = 0.5*(interval_default_max_-interval_default_min_);
 interval_range_ = interval_max_-interval_min_;
 interval_mean_  = 0.5*(interval_max_-interval_min_);
 argT_derivf_ = interval_default_range_/interval_range_;
}

double BasisFunctions::translateArgument(const double arg, bool inside_range)
{
 inside_range=true;
 double argT = (arg-interval_mean_)*argT_derivf_;
 if(argT < interval_default_min_){
  inside_range=false;
  argT=interval_default_min_;
 }
 else if(argT > interval_default_max_){
  inside_range=false;
  argT=interval_default_max_;
 }
 return argT;
}

void BasisFunctions::apply(){}

void BasisFunctions::calculate(){}

void BasisFunctions::setupBF()
{
 if(interval_default_min_>interval_default_max_){error("setupBF: default intervals are not correctly set");}
 setupInterval();
 setupDescription();
 if(bf_description_.size()==1){error("setupBF: the description of the basis functions is not correct.");}
 setupBFIntegrals();
 if(bf_integrals_.size()==1){error("setupBF: the integrals of the basis functions is not correct.");}
 if(type_=="Undefined"){error("setupBF: the type of the basis function is not defined.");}
 has_been_set=true;
}

void BasisFunctions::printInfo()
{
 if(!has_been_set){error("the basis set has not be setup correctly");}
 log.printf("  One-dimensional basis set\n");
 log.printf("   Type: %s\n",type_.c_str());
 if(periodic_){log.printf("   The basis functions are periodic\n");}
 log.printf("   Default interval of basis set: %f to %f\n",interval_default_min_,interval_default_max_);
 log.printf("   Defined interval of basis set: %f to %f\n",interval_min_,interval_max_);
 log.printf("   Derivative factor due to interval translation: %f\n",argT_derivf_);
 log.printf("   Order of basis set: %d\n",norder_);
 log.printf("   Number of basis functions: %d\n",nbasis_);
 log.printf("   Basis functions:\n");
 for(unsigned int i=0; i < nbasis_;i++){log.printf("    %2d       %10s\n",i,bf_description_[i].c_str());}
 log.printf("   --------------------------\n");
 log.printf("   Basis functions integrals:\n");
 for(unsigned int i=0; i < nbasis_;i++){log.printf("    %2d       %f\n",i,bf_integrals_[i]);}
 log.printf("   --------------------------\n");
}

}
}


