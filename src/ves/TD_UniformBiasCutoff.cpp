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
#include "tools/Grid.h"


namespace PLMD{
namespace ves{

//+PLUMEDOC VES_TARGETDIST_HIDDEN UNIFORM_BIAS_CUTOFF
/*
Uniform target distribution with a bias cutoff (dynamic).

\par Examples

*/
//+ENDPLUMEDOC


class TD_UniformBiasCutoff : public TargetDistribution {
private:
public:
  static void registerKeywords( Keywords&);
  explicit TD_UniformBiasCutoff( const TargetDistributionOptions& to );
  void updateGrid();
  double getValue(const std::vector<double>&) const;
};


VES_REGISTER_TARGET_DISTRIBUTION(TD_UniformBiasCutoff,"UNIFORM_BIAS_CUTOFF")


void TD_UniformBiasCutoff::registerKeywords(Keywords& keys) {
  TargetDistribution::registerKeywords(keys);
  keys.use("BIAS_CUTOFF");
}


TD_UniformBiasCutoff::TD_UniformBiasCutoff(const TargetDistributionOptions& to):
TargetDistribution(to)
{
  if(!biasCutoffActive()){
    plumed_merror("using UNIFORM_BIAS_CUTOFF without a BIAS_CUTOFF keywords does not make sense");
  }
  checkRead();
}


double TD_UniformBiasCutoff::getValue(const std::vector<double>& argument) const {
  plumed_merror("getValue not implemented for TD_UniformBiasCutoff");
  return 0.0;
}


void TD_UniformBiasCutoff::updateGrid(){
  plumed_massert(biasCutoffActive(),"The UNIFORM_BIAS_CUTOFF should only be used when employing a bias cutoff");
  for(Grid::index_t l=0; l<targetDistGrid().getSize(); l++){
    targetDistGrid().setValue(l,1.0);
  }
  logTargetDistGrid().clear();
}


}
}
