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
  double inverse_normalization;
public:
  static void registerKeywords( Keywords&);
  UniformDistribution( const TargetDistributionOptions& to );
  double distribution(const std::vector<double> argument);
  double getNormalization() const;
};


VARIATIONAL_REGISTER_TARGET_DISTRIBUTION(UniformDistribution,"UNIFORM")


void UniformDistribution::registerKeywords(Keywords& keys) {
  TargetDistributionBase::registerKeywords(keys);
  keys.add("optional","NORMALIZATION","Normalization factor for the uniform distribution. If not given the distribution will be left unnormalized.");
}


UniformDistribution::UniformDistribution(const TargetDistributionOptions& to):
TargetDistributionBase(to)
{
  double normalization;
  if(parse("NORMALIZATION",normalization,true)){
    inverse_normalization=1.0/normalization;
    setNormalized();
  }
  else{
    inverse_normalization=1.0;
    setNotNormalized();
  }
  checkRead();
}


double UniformDistribution::distribution(const std::vector<double> argument) {
  double value=inverse_normalization;
  return value;
}


double UniformDistribution::getNormalization() const {
  return 1.0/inverse_normalization;
}



}
