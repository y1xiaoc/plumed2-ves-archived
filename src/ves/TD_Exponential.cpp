/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2016-2017 The ves-code team
   (see the PEOPLE-VES file at the root of this folder for a list of names)

   See http://www.ves-code.org for more information.

   This file is part of ves-code, version 1.

   ves-code is free software: you can redistribute it and/or modify
   it under the terms of the GNU Lesser General Public License as published by
   the Free Software Foundation, either version 3 of the License, or
   (at your option) any later version.

   ves-code is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU Lesser General Public License for more details.

   You should have received a copy of the GNU Lesser General Public License
   along with ves-code.  If not, see <http://www.gnu.org/licenses/>.
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */

#include "TargetDistribution.h"
#include "TargetDistributionRegister.h"

#include "tools/Keywords.h"


namespace PLMD{
namespace ves{

//+PLUMEDOC VES_TARGETDIST EXPONENTIAL HIDDEN
/*
Exponential distribution (static).

\par Examples

*/
//+ENDPLUMEDOC

class TD_Exponential: public TargetDistribution {
  std::vector<double> minima_;
  std::vector<double> lambda_;
public:
  static void registerKeywords(Keywords&);
  explicit TD_Exponential(const TargetDistributionOptions& to);
  double getValue(const std::vector<double>&) const;
};


VES_REGISTER_TARGET_DISTRIBUTION(TD_Exponential,"EXPONENTIAL")


void TD_Exponential::registerKeywords(Keywords& keys){
  TargetDistribution::registerKeywords(keys);
  keys.add("compulsory","MINIMA","The minima of the exponential distribution.");
  keys.add("compulsory","LAMBDA","The lambda parameters for the chi-squared distribution.");
}


TD_Exponential::TD_Exponential( const TargetDistributionOptions& to ):
TargetDistribution(to),
minima_(0),
lambda_(0)
{
  parseVector("MINIMA",minima_);
  parseVector("LAMBDA",lambda_);
  for(unsigned int k=0; k<lambda_.size(); k++){
    if(lambda_[k] < 0.0){plumed_merror(getName()+": the values given in LAMBDA should be postive.");}
  }


  setDimension(minima_.size());
  if(lambda_.size()!=getDimension()){plumed_merror(getName()+": the LAMBDA keyword does not match the given dimension in MINIMA");}
  checkRead();
}


double TD_Exponential::getValue(const std::vector<double>& argument) const {
  double value = 1.0;
  for(unsigned int k=0; k<argument.size(); k++){
    double arg = (argument[k]-minima_[k])*lambda_[k];
    if(arg<0.0){plumed_merror(getName()+": the exponential distribution is not defined for values less that ones given in MINIMA");}
    value *= lambda_[k]*exp(-arg);
  }
  return value;
}



}
}
