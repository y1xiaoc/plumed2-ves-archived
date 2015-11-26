/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2011-2015 The plumed team
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
#include "VesBias.h"
#include "ves_tools/CoeffsVector.h"
#include "ves_tools/CoeffsMatrix.h"
#include "ves_basisfunctions/BasisFunctions.h"

#include "bias/Bias.h"
#include "core/ActionRegister.h"
#include "core/ActionSet.h"
#include "core/PlumedMain.h"


using namespace std;


namespace PLMD{
namespace bias{

//+PLUMEDOC BIAS MOVINGRESTRAINT
/*

*/
//+ENDPLUMEDOC


class TestVesBias : public VesBias{
private:
  BasisFunctions* bf_pointer;
  Value* valueBias;
  Value* valueForce2;
public:
  explicit TestVesBias(const ActionOptions&);
  void calculate();
  static void registerKeywords( Keywords& keys );
};

PLUMED_REGISTER_ACTION(TestVesBias,"TEST_VES_BIAS")

void TestVesBias::registerKeywords( Keywords& keys ){
  VesBias::registerKeywords(keys);
  keys.use("ARG");
  keys.add("compulsory","BASIS_SET","the label of the basis set that you want to use");
}

TestVesBias::TestVesBias(const ActionOptions&ao):
PLUMED_VESBIAS_INIT(ao)
{
  std::string basisset_label="";
  parse("BASIS_SET",basisset_label);
  checkRead();
  bf_pointer=plumed.getActionSet().selectWithLabel<BasisFunctions*>(basisset_label);
  std::vector<BasisFunctions*> bf(1);
  bf[0]=bf_pointer;
  std::vector<Value*> args(1);
  args[0]=getArguments()[0];
  initializeCoeffs(args,bf);
  setCoeffsDerivsOverTargetDist(bf_pointer->getBasisFunctionIntegrals());
  addComponent("bias"); componentIsNotPeriodic("bias");
  addComponent("force2"); componentIsNotPeriodic("force2");
  valueBias=getPntrToComponent("bias");
  valueForce2=getPntrToComponent("force2");
}


void TestVesBias::calculate() {
  std::vector<double> bf_values(bf_pointer->getSize());
  std::vector<double> bf_derivs(bf_pointer->getSize());
  bool inside=true;
  double cv = getArgument(0);
  double cvT=0.0;
  bf_pointer->getAllValues(cv,cvT,inside,bf_values,bf_derivs);
  double bias = 0.0;
  double deriv = 0.0;
  for(size_t i=0; i<numberOfCoeffs(); i++){
    bias += Coeffs()[i]*bf_values[i];
    deriv += -1.0*Coeffs()[i]*bf_derivs[i];
  }
  valueBias->set(bias);
  valueForce2->set(deriv*deriv);
  setOutputForce(0,deriv);
  setCoeffsDerivs(bf_values);
}

}
}
