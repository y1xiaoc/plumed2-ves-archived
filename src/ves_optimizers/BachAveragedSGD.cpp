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
#include "ves_tools/CoeffsMatrix.h"

#include "core/ActionRegister.h"


namespace PLMD{

class BachAveragedSGD : public Optimizer {
private:
public:
  static void registerKeywords(Keywords&);
  explicit BachAveragedSGD(const ActionOptions&);
  void coeffsUpdate();
};


PLUMED_REGISTER_ACTION(BachAveragedSGD,"AVERAGED_SGD")


void BachAveragedSGD::registerKeywords(Keywords& keys){
  Optimizer::registerKeywords(keys);
  //
  Optimizer::useFixedStepSizeKeywords(keys);
  Optimizer::useHessianKeywords(keys);
  keys.use("MASK_FILE");
}


BachAveragedSGD::BachAveragedSGD(const ActionOptions&ao):
PLUMED_OPTIMIZER_INIT(ao)
{
  turnOnHessian();
  checkRead();
}


void BachAveragedSGD::coeffsUpdate() {
  double aver_decay = 1.0 / ( getIterationCounterDbl() + 1.0 );
  AuxCoeffs() = AuxCoeffs() - StepSize()*CoeffsMask() * ( Gradient() + Hessian()*(AuxCoeffs()-Coeffs()) );
  //AuxCoeffs() = AuxCoeffs() - StepSize() * ( Gradient() + Hessian()*(AuxCoeffs()-Coeffs()) );
  Coeffs() += aver_decay * ( AuxCoeffs()-Coeffs() );
}


}