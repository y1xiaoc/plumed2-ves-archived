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
#include "TargetDistribution.h"
#include "TargetDistributionRegister.h"

#include "tools/Keywords.h"

namespace PLMD {

class UniformDistribution : public TargetDistribution {
  double normalization_;
  double inverse_normalization_;
  std::vector<double> minima_;
  std::vector<double> maxima_;
public:
  static void registerKeywords( Keywords&);
  explicit UniformDistribution( const TargetDistributionOptions& to );
  double getValue(const std::vector<double>&) const;
  double getNormalization() const {return normalization_;}
};


VARIATIONAL_REGISTER_TARGET_DISTRIBUTION(UniformDistribution,"UNIFORM")


void UniformDistribution::registerKeywords(Keywords& keys) {
  TargetDistribution::registerKeywords(keys);
  keys.add("compulsory","MINIMA","The minima of the intervals on which the target distribution is defined.");
  keys.add("compulsory","MAXIMA","The maxima of the intervals on which the target distribution is defined.");
}


UniformDistribution::UniformDistribution(const TargetDistributionOptions& to):
TargetDistribution(to),
normalization_(1.0),
inverse_normalization_(1.0),
minima_(0),
maxima_(0)
{
  parseVector("MINIMA",minima_);
  parseVector("MAXIMA",maxima_);
  plumed_massert(minima_.size()==maxima_.size(),"MINIMA and MAXIMA for the uniform distribution do not have the same size");
  setDimension(minima_.size());
  normalization_ = 1.0;
  for(unsigned int k=0; k<getDimension(); k++){
    plumed_massert(maxima_[k]>minima_[k],"Check MINIMA and MAXIMA keywords");
    normalization_ *= maxima_[k]-minima_[k];
  }
  inverse_normalization_=1.0/normalization_;
  setNormalized();
  checkRead();
}


double UniformDistribution::getValue(const std::vector<double>& argument) const {
  double outside = 0.0;
  double inside = inverse_normalization_;
  for(unsigned int k=0; k<getDimension(); k++){
    if(argument[k] < minima_[k] || argument[k] > maxima_[k]){return outside;}
  }
  return inside;
}


}
