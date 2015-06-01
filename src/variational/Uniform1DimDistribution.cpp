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
#include "TargetDistribution1DimBase.h"
#include "TargetDistribution1DimRegister.h"
#include "tools/Keywords.h"

namespace PLMD {

class Uniform1DimDistribution : public TargetDistribution1DimBase {
 double inverse_normalization;
public:
  static void registerKeywords( Keywords&);
  Uniform1DimDistribution( const TargetDistribution1DimOptions& to );
  double distribution(const double);
};

VARIATIONAL_REGISTER_TARGET_DISTRIBUTION_1D(Uniform1DimDistribution,"UNIFORM")

void Uniform1DimDistribution::registerKeywords( Keywords& keys )
{
 TargetDistribution1DimBase::registerKeywords(keys);
 keys.add("optional","NORMALIZATION","Normalization factor for the uniform distribution. If not given the distribution will be left unnormalized.");
}

Uniform1DimDistribution::Uniform1DimDistribution( const TargetDistribution1DimOptions& to ):
TargetDistribution1DimBase(to)
{
 double normalization;
 if(parse("NORMALIZATION",normalization,true)){
  inverse_normalization=1.0/normalization;
  setNormalized();
 }else{
  inverse_normalization=1.0;
  setNotNormalized();
 }
}

double Uniform1DimDistribution::distribution(const double argument)
{
 double value=inverse_normalization;
 return value;
}

}
