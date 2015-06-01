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
#include "BasisFunctions.h"

namespace PLMD {
namespace variational {

class Uniform1DimDistribution : public TargetDistribution1DimBase {
public:
  Uniform1DimDistribution( const TargetDistribution1DimOptions& to );
  double calculate_ps(double);
};

VARIATIONAL_REGISTER_TARGET_DISTRIBUTION_1D(Uniform1DimDistribution,"UNIFORM")

Uniform1DimDistribution::Uniform1DimDistribution( const TargetDistribution1DimOptions& to ):
TargetDistribution1DimBase(to)
{
 setNormalizationFactor(BF_pointer()->intervalRange());
}

double Uniform1DimDistribution::calculate_ps(double cv_value)
{
 double ps_value=1.0/getNormalizationFactor();
 return ps_value;
}

}
}
