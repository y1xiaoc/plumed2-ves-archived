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
#include "TargetDistributionBase.h"
#include "TargetDistributionRegister.h"

#include "tools/Keywords.h"

namespace PLMD {

class UniformDistribution : public TargetDistributionBase {
  double normalization;
  double inverse_normalization;
  std::vector<double> minima;
  std::vector<double> maxima;
public:
  static void registerKeywords( Keywords&);
  explicit UniformDistribution( const TargetDistributionOptions& to );
  double getValue(const std::vector<double>&) const;
  double getNormalization() const;
};


VARIATIONAL_REGISTER_TARGET_DISTRIBUTION(UniformDistribution,"UNIFORM")


void UniformDistribution::registerKeywords(Keywords& keys) {
  TargetDistributionBase::registerKeywords(keys);
  keys.add("compulsory","MINIMA","The minima of the intervals on which the target distribution is defined.");
  keys.add("compulsory","MAXIMA","The maxima of the intervals on which the target distribution is defined.");
}


UniformDistribution::UniformDistribution(const TargetDistributionOptions& to):
TargetDistributionBase(to)
{
  parseVector("MINIMA",minima);
  parseVector("MAXIMA",maxima);
  plumed_massert(minima.size()==maxima.size(),"MINIMA and MAXIMA for the uniform distribution do not have the same size");
  setDimension(minima.size());
  normalization = 1.0;
  for(unsigned int k=0; k<getDimension(); k++){
    plumed_massert(maxima[k]>minima[k],"Check MINIMA and MAXIMA keywords");
    normalization *= maxima[k]-minima[k];
  }
  inverse_normalization=1.0/normalization;
  setNormalized();
  checkRead();
}


double UniformDistribution::getValue(const std::vector<double>& argument) const {
  double outside = 0.0;
  double inside = inverse_normalization;
  for(unsigned int k=0; k<getDimension(); k++){
    if(argument[k] < minima[k] || argument[k] > maxima[k]){return outside;}
  }
  return inside;
}


double UniformDistribution::getNormalization() const {
  return normalization;
}



}
