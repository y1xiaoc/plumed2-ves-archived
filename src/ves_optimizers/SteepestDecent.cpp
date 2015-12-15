/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2011-2014 The plumed team
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
#include "Optimizer.h"
#include "ves_tools/CoeffsVector.h"

#include "core/ActionRegister.h"


namespace PLMD{

class SteepestDecent : public Optimizer {

public:
  static void registerKeywords(Keywords&);
  explicit SteepestDecent(const ActionOptions&);
  void coeffsUpdate();
};


PLUMED_REGISTER_ACTION(SteepestDecent,"STEEPEST_DECENT")


void SteepestDecent::registerKeywords(Keywords& keys){
  Optimizer::registerKeywords(keys);
  Optimizer::useFixedStepSizeKeywords(keys);
  Optimizer::useMultipleWalkersKeywords(keys);
}


SteepestDecent::SteepestDecent(const ActionOptions&ao):
PLUMED_OPTIMIZER_INIT(ao)
{
  checkRead();
}


void SteepestDecent::coeffsUpdate() {
  Coeffs() = Coeffs() - StepSize()*Gradient();
}


}
