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
#include "tools/Grid.h"

namespace PLMD {

class UniformBiasCutoffDistribution : public TargetDistribution {
public:
  static void registerKeywords( Keywords&);
  explicit UniformBiasCutoffDistribution( const TargetDistributionOptions& to );
  void updateGrid();
  double getValue(const std::vector<double>&) const;
};


VES_REGISTER_TARGET_DISTRIBUTION(UniformBiasCutoffDistribution,"UNIFORM_BIAS_CUTOFF")


void UniformBiasCutoffDistribution::registerKeywords(Keywords& keys) {
  TargetDistribution::registerKeywords(keys);
}


UniformBiasCutoffDistribution::UniformBiasCutoffDistribution(const TargetDistributionOptions& to):
TargetDistribution(to)
{
  setNormalized();
  setupBiasCutoff();
  checkRead();
}


double UniformBiasCutoffDistribution::getValue(const std::vector<double>& argument) const {
  plumed_merror("getValue not implemented for UniformBiasCutoffDistribution");
  return 0.0;
}


void UniformBiasCutoffDistribution::updateGrid(){
  plumed_massert(biasCutoffActive(),"The UNIFORM_BIAS_CUTOFF should only be used when employing a bias cutoff");
  for(Grid::index_t l=0; l<targetDistGrid().getSize(); l++){
    targetDistGrid().setValue(l,1.0);
  }
}


}
