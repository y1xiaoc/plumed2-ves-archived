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


namespace PLMD {
namespace ves {

//+PLUMEDOC VES_TARGETDIST_HIDDEN CHI_SQUARED
/*
Chi-squared distribution (static).

\par Examples

*/
//+ENDPLUMEDOC

class TD_ChiSquared: public TargetDistribution {
  std::vector<double> minima_;
  std::vector<double> sigma_;
  std::vector<double> kappa_;
  std::vector<double> normalization_;
public:
  static void registerKeywords(Keywords&);
  explicit TD_ChiSquared(const TargetDistributionOptions& to);
  double getValue(const std::vector<double>&) const;
};


VES_REGISTER_TARGET_DISTRIBUTION(TD_ChiSquared,"CHI_SQUARED")


void TD_ChiSquared::registerKeywords(Keywords& keys) {
  TargetDistribution::registerKeywords(keys);
  keys.add("compulsory","MINIMA","The minima of the chi-squared distribution.");
  keys.add("compulsory","SIGMA","The sigma parameters for the chi-squared distribution.");
  keys.add("compulsory","KAPPA","The kappa parameters for the chi-squared distribution.");
  keys.use("BIAS_CUTOFF");
  keys.use("WELLTEMPERED_FACTOR");
  keys.use("SHIFT_TO_ZERO");
  keys.use("NORMALIZE");
}


TD_ChiSquared::TD_ChiSquared( const TargetDistributionOptions& to ):
  TargetDistribution(to),
  minima_(0),
  sigma_(0),
  kappa_(0),
  normalization_(0)
{
  parseVector("MINIMA",minima_);
  parseVector("SIGMA",sigma_);
  for(unsigned int k=0; k<sigma_.size(); k++) {
    if(sigma_[k] < 0.0) {plumed_merror(getName()+": the values given in SIGMA should be postive.");}
  }

  std::vector<unsigned int> kappa_int(0);
  parseVector("KAPPA",kappa_int,true);
  if(kappa_int.size()==0) {plumed_merror(getName()+": some problem with KAPPA keyword, should given as postive integer(s) larger than 1");}
  kappa_.resize(kappa_int.size());
  for(unsigned int k=0; k<kappa_int.size(); k++) {
    if(kappa_int[k] < 2) {plumed_merror(getName()+": KAPPA should be a integers 2 or higher");}
    kappa_[k] = static_cast<double>(kappa_int[k]);
  }

  setDimension(minima_.size());
  if(sigma_.size()!=getDimension()) {plumed_merror(getName()+": the SIGMA keyword does not match the given dimension in MINIMA");}
  if(kappa_.size()!=getDimension()) {plumed_merror(getName()+": the KAPPA keyword does not match the given dimension in MINIMA");}

  normalization_.resize(getDimension());
  for(unsigned int k=0; k<getDimension(); k++) {
    normalization_[k] = 1.0/(pow(2.0,0.5*kappa_[k])*tgamma(0.5*kappa_[k])*sigma_[k]);
  }
  checkRead();
}


double TD_ChiSquared::getValue(const std::vector<double>& argument) const {
  double value = 1.0;
  for(unsigned int k=0; k<argument.size(); k++) {
    double arg=(argument[k]-minima_[k])/sigma_[k];
    if(arg<0.0) {plumed_merror(getName()+": the chi-squared istribution is not defined for values less that ones given in MINIMA");}
    value *= normalization_[k] * pow(arg,0.5*kappa_[k]-1.0) * exp(-0.5*arg);
  }
  return value;
}


}
}
