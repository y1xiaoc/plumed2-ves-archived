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

//+PLUMEDOC INTERNAL UNIFORM
/*
  Uniform target distribution
*/
//+ENDPLUMEDOC

class TD_Uniform : public TargetDistribution {
  double normalization_;
  double inverse_normalization_;
  std::vector<double> minima_;
  std::vector<double> maxima_;
public:
  static void registerKeywords( Keywords&);
  explicit TD_Uniform( const TargetDistributionOptions& to );
  double getValue(const std::vector<double>&) const;
};


VES_REGISTER_TARGET_DISTRIBUTION(TD_Uniform,"UNIFORM")


void TD_Uniform::registerKeywords(Keywords& keys) {
  TargetDistribution::registerKeywords(keys);
  keys.add("compulsory","MINIMA","The minima of the intervals on which the target distribution is defined.");
  keys.add("compulsory","MAXIMA","The maxima of the intervals on which the target distribution is defined.");
}


TD_Uniform::TD_Uniform(const TargetDistributionOptions& to):
TargetDistribution(to),
normalization_(1.0),
inverse_normalization_(1.0),
minima_(0),
maxima_(0)
{
  std::vector<std::string> min_str_;
  std::vector<std::string> max_str_;
  parseVector("MINIMA",min_str_);
  parseVector("MAXIMA",max_str_);
  plumed_massert(min_str_.size()==max_str_.size(),"MINIMA and MAXIMA for the uniform distribution do not have the same size");
  //
  setDimension(min_str_.size());
  minima_.assign(getDimension(),0.0);
  maxima_.assign(getDimension(),0.0);
  //
  normalization_ = 1.0;
  for(unsigned int k=0; k<getDimension(); k++){
    if(!Tools::convert(min_str_[k],minima_[k])){
      plumed_merror("cannot convert one of the values given in MINIMA to a double");
    }
    if(!Tools::convert(max_str_[k],maxima_[k])){
      plumed_merror("cannot convert one of the values given in MAXIMA to a double");
    }
    plumed_massert(maxima_[k]>minima_[k],"Check MINIMA and MAXIMA keywords");
    normalization_ *= maxima_[k]-minima_[k];
  }
  inverse_normalization_=1.0/normalization_;
  checkRead();
}


double TD_Uniform::getValue(const std::vector<double>& argument) const {
  double outside = 0.0;
  double inside = inverse_normalization_;
  for(unsigned int k=0; k<getDimension(); k++){
    if(argument[k] < minima_[k] || argument[k] > maxima_[k]){return outside;}
  }
  return inside;
}


}
}
