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
  std::vector<BasisFunctions*> bf_pointers;
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
  keys.add("compulsory","BASIS_SET","the label of the basis sets that you want to use");
}

TestVesBias::TestVesBias(const ActionOptions&ao):
PLUMED_VESBIAS_INIT(ao),
bf_pointers(getNumberOfArguments(),NULL)
{
  std::vector<std::string> basisset_labels;
  parseVector("BASIS_SET",basisset_labels);
  plumed_massert(basisset_labels.size()==getNumberOfArguments(),"number of arguments should match the number of basis set labels");
  checkRead();

  for(unsigned int i=0; i<basisset_labels.size(); i++){
    bf_pointers[i] = plumed.getActionSet().selectWithLabel<BasisFunctions*>(basisset_labels[i]);
    plumed_massert(bf_pointers[i]!=NULL,"error in basis set");
  }

  if(getNumberOfArguments()>1){enableMultipleCoeffsSets();}

  for(unsigned int i=0; i<getNumberOfArguments(); i++){
    std::vector<Value*> arg(1);
    arg[0]=getArguments()[i];
    std::vector<BasisFunctions*> bf(1);
    bf[0]=bf_pointers[i];
    addCoeffsSet(arg,bf);
  }

  for(unsigned int i=0; i<numberOfCoeffsSets(); i++){
    setCoeffsDerivsOverTargetDist(bf_pointers[i]->getBasisFunctionIntegrals(),i);
  }

  addComponent("bias"); componentIsNotPeriodic("bias");
  valueBias=getPntrToComponent("bias");

}


void TestVesBias::calculate() {
  double total_bias = 0.0;
  for(unsigned int k=0; k<getNumberOfArguments(); k++){
    std::vector<double> bf_values(bf_pointers[k]->getSize());
    std::vector<double> bf_derivs(bf_pointers[k]->getSize());
    bool inside=true;
    double cv = getArgument(k);
    double cvT=0.0;
    bf_pointers[k]->getAllValues(cv,cvT,inside,bf_values,bf_derivs);
    double bias = 0.0;
    double deriv = 0.0;
    for(size_t i=0; i<numberOfCoeffs(k); i++){
      bias += Coeffs(k)[i]*bf_values[i];
      deriv += -1.0*Coeffs(k)[i]*bf_derivs[i];
    }
    total_bias += bias;
    setOutputForce(k,deriv);
    setCoeffsDerivs(bf_values,k);
  }
  valueBias->set(total_bias);
}

}
}
