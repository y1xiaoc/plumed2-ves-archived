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

//+PLUMEDOC VES_TARGETDIST UNIFORM
/*
Uniform target distribution (static).

\par Examples

*/
//+ENDPLUMEDOC

class TD_UniformWithSwitchingFunction : public TargetDistribution {
  std::vector<double> minima_;
  std::vector<double> maxima_;
  std::vector<double> sigma_min_;
  std::vector<double> sigma_max_;
  double GaussianSwitchingFunc(const double, const double, const double) const;
public:
  static void registerKeywords( Keywords&);
  explicit TD_UniformWithSwitchingFunction( const TargetDistributionOptions& to );
  double getValue(const std::vector<double>&) const;
};


VES_REGISTER_TARGET_DISTRIBUTION(TD_UniformWithSwitchingFunction,"UNIFORM")


void TD_UniformWithSwitchingFunction::registerKeywords(Keywords& keys) {
  TargetDistribution::registerKeywords(keys);
  keys.add("optional","MINIMA","The minima of the intervals where the target distribution is uniform.");
  keys.add("optional","MAXIMA","The maxima of the intervals where the target distribution is uniform.");
  keys.add("optional","SIGMA_MINIMA","The sigma values of the Gaussian switching functions for the minima of the intervals. Value of 0.0 means that switch is done without a smooth switching function, this is the default behaviour.");
  keys.add("optional","SIGMA_MAXIMA","The sigma values of the Gaussian switching functions for the maxima of the intervals. Value of 0.0 means that switch is done without a smooth switching function, this is the default behaviour.");
  keys.use("BIAS_CUTOFF");
}


TD_UniformWithSwitchingFunction::TD_UniformWithSwitchingFunction(const TargetDistributionOptions& to):
TargetDistribution(to),
minima_(0),
maxima_(0),
sigma_min_(0),
sigma_max_(0)
{
  parseVector("MINIMA",minima_,true);
  parseVector("MAXIMA",maxima_,true);
  if(minima_.size()!=maxima_.size()){
    plumed_merror(getName()+": MINIMA and MAXIMA do not have the same size");
  }
  setDimension(minima_.size());
  for(unsigned int k=0; k<getDimension(); k++){
    if(minima_[k]>maxima_[k]){
      plumed_merror(getName()+": error in MINIMA and MAXIMA keywords, one of the MINIMA values is larger than the corresponding MAXIMA values");
    }
  }
  //
  parseVector("SIGMA_MINIMA",sigma_min_,true);
  parseVector("SIGMA_MAXIMA",sigma_max_,true);
  if(sigma_min_.size()==0){sigma_min_.assign(getDimension(),0.0);}
  if(sigma_max_.size()==0){sigma_max_.assign(getDimension(),0.0);}
  if(sigma_min_.size()!=getDimension()){plumed_merror(getName()+": SIGMA_MINIMA has the wrong size");}
  if(sigma_max_.size()!=getDimension()){plumed_merror(getName()+": SIGMA_MAXIMA has the wrong size");}
  //
  setForcedNormalization();
  checkRead();
}


double TD_UniformWithSwitchingFunction::getValue(const std::vector<double>& argument) const {
  //
  if(minima_.size()==0){return 1.0;}
  //
  double value = 1.0;
  for(unsigned int k=0; k<getDimension(); k++){
    double tmp;
    if(argument[k] < minima_[k]){
      tmp = GaussianSwitchingFunc(argument[k],minima_[k],sigma_min_[k]);
    }
    else if(argument[k] > maxima_[k]){
      tmp = GaussianSwitchingFunc(argument[k],maxima_[k],sigma_max_[k]);
    }
    else {
      tmp = 1.0;
    }
    value *= tmp;
  }
  return value;
}

inline
double TD_UniformWithSwitchingFunction::GaussianSwitchingFunc(const double argument, const double center, const double sigma) const {
  if(sigma>0.0){
    double arg=(argument-center)/sigma;
    return exp(-0.5*arg*arg);
  }
  else{
    return 0.0;
  }
}






}
}
