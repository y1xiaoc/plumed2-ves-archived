/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2015-2016 The ves-code team
   (see the PEOPLE-VES file at the root of the distribution for a list of names)

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

//+PLUMEDOC VES_TARGETDIST GENERALIZED_EXTREME_VALUE
/*
Generalized extreme value distribution (static).

\par Examples

*/
//+ENDPLUMEDOC

class TD_GeneralizedExtremeValue: public TargetDistribution {
  std::vector<double> center_;
  std::vector<double> sigma_;
  std::vector<double> epsilon_;
  std::vector<double> normalization_;
  double GEVdiagonal(const std::vector<double>&, const std::vector<double>&, const std::vector<double>&, const std::vector<double>&, const std::vector<double>&) const;
public:
  static void registerKeywords(Keywords&);
  explicit TD_GeneralizedExtremeValue(const TargetDistributionOptions& to);
  double getValue(const std::vector<double>&) const;
};


VES_REGISTER_TARGET_DISTRIBUTION(TD_GeneralizedExtremeValue,"GENERALIZED_EXTREME_VALUE")


void TD_GeneralizedExtremeValue::registerKeywords(Keywords& keys){
  TargetDistribution::registerKeywords(keys);
  keys.add("compulsory","CENTER","The center of the generalized extreme value distribution.");
  keys.add("compulsory","SIGMA","The sigma (scale) parameters for the generalized extreme value distribution.");
  keys.add("compulsory","EPSILON","The epsilon (shape) parameters for the generalized extreme value distribution.");
}


TD_GeneralizedExtremeValue::TD_GeneralizedExtremeValue( const TargetDistributionOptions& to ):
TargetDistribution(to),
center_(0),
sigma_(0),
epsilon_(0),
normalization_(0)
{
  parseVector("CENTER",center_);
  parseVector("SIGMA",sigma_);
  parseVector("EPSILON",epsilon_);

  setDimension(center_.size());
  if(sigma_.size()!=getDimension()){plumed_merror(getName()+": the SIGMA keyword does not match the given dimension in MINIMA");}
  if(epsilon_.size()!=getDimension()){plumed_merror(getName()+": the EPSILON keyword does not match the given dimension in MINIMA");}

  normalization_.resize(getDimension());
  for(unsigned int k=0; k<getDimension(); k++){
    if(sigma_[k]<0.0){plumed_merror(getName()+": the values given in SIGMA should be larger then 0.0");}
    normalization_[k] = 1.0/sigma_[k];
  }
  checkRead();
}


double TD_GeneralizedExtremeValue::getValue(const std::vector<double>& argument) const {
  return GEVdiagonal(argument,center_,sigma_,epsilon_,normalization_);
}


double TD_GeneralizedExtremeValue::GEVdiagonal(const std::vector<double>& argument, const std::vector<double>& center, const std::vector<double>& sigma, const std::vector<double>& epsilon, const std::vector<double>& normalization) const {
  double value = 1.0;
  for(unsigned int k=0; k<argument.size(); k++){
    double arg=(argument[k]-center[k])/sigma[k];
    double tx;
    if(epsilon_[k]!=0.0){
      if( epsilon_[k]>0 && argument[k] <= (center[k]-sigma[k]/epsilon[k]) ){return 0.0;}
      if( epsilon_[k]<0 && argument[k] > (center[k]-sigma[k]/epsilon[k]) ){return 0.0;}
      tx = pow( (1.0+arg*epsilon[k]) , -1.0/epsilon[k] );
    }
    else{
      tx = exp(-arg);
    }
    value *= normalization[k] * pow(tx,epsilon[k]+1.0) * exp(-tx);
  }
  return value;
}



}
}
