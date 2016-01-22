/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2014 The plumed team
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
#include "ves_tools/CoeffsVector.h"
#include "ves_tools/CoeffsMatrix.h"
#include "ves_targetdistributions/TargetDistributionBase.h"
#include "BasisFunctions.h"

#include "../function/Function.h"
#include "core/ActionRegister.h"
#include "core/ActionSet.h"
#include "core/PlumedMain.h"
#include "tools/File.h"
#include "tools/Communicator.h"



// using namespace std;

namespace PLMD{
namespace function{

//+PLUMEDOC FUNCTION TEST_BASISFUNCTIONS
/*

*/
//+ENDPLUMEDOC


class TestBFs :
  public Function
{
  BasisFunctions* bf_pointer;
  BasisFunctions* bf_pointer2;
  CoeffsVector* coeffs2;
  unsigned int bf_order_;
public:
  explicit TestBFs(const ActionOptions&);
  void calculate();
  static void registerKeywords(Keywords& keys);
};


PLUMED_REGISTER_ACTION(TestBFs,"TEST_BASISFUNCTIONS")

void TestBFs::registerKeywords(Keywords& keys){
  Function::registerKeywords(keys);
  keys.addOutputComponent("value","default","value of basis function");
  keys.addOutputComponent("deriv","default","deriative of basis function");
  keys.addOutputComponent("arg","default","input argument");
  keys.addOutputComponent("argT","default","translated input argument");
  keys.addOutputComponent("inside","default","1.0 if inside interval, otherwise 0.0");
  keys.use("ARG");
  keys.add("compulsory","BASIS_SET","the label of the basis set that you want to use");
  keys.add("compulsory","BASIS_SET2","the label of the basis set that you want to use");
  keys.add("compulsory","N","the number of the basis function you want to test");
}

TestBFs::TestBFs(const ActionOptions&ao):
Action(ao),
Function(ao)
{
  // if(getNumberOfArguments()>1){error("only one argument allowed");}
  std::string basisset_label="";
  parse("BASIS_SET",basisset_label);
  bf_pointer=plumed.getActionSet().selectWithLabel<BasisFunctions*>(basisset_label);
  parse("BASIS_SET2",basisset_label);
  bf_pointer2=plumed.getActionSet().selectWithLabel<BasisFunctions*>(basisset_label);
  // bf_pointer->printInfo();
  parse("N",bf_order_);
  addComponent("value"); componentIsNotPeriodic("value");
  addComponent("deriv"); componentIsNotPeriodic("deriv");
  addComponent("arg"); componentIsNotPeriodic("arg");
  addComponent("argT"); componentIsNotPeriodic("argT");
  addComponent("inside"); componentIsNotPeriodic("inside");
  checkRead();
  log.printf("  using the %u order basis function from the %s basis set\n",bf_order_,basisset_label.c_str());

  std::vector<BasisFunctions*> bf; bf.resize(2); bf[0]=bf_pointer; bf[1]=bf_pointer2;
  std::vector<BasisFunctions*> bf2; bf2.resize(2); bf2[0]=bf_pointer2; bf2[1]=bf_pointer;
  std::vector<Value*> args; args.resize(2); args[0]=getArguments()[0]; args[1]=getArguments()[1];

  CoeffsMatrix coeffsM = CoeffsMatrix("coeffsM",args,bf,comm);
  coeffsM.randomizeValuesGaussian(1);
  coeffsM.writeToFile("coeffsM.data");

  CoeffsVector coeffsM2V = CoeffsVector("coeffsM2V",&coeffsM,comm);
  coeffsM2V = 3.0; coeffsM2V[4]=414.0;
  coeffsM2V.writeToFile("coeffsM2V.data");


  std::vector<std::vector<BasisFunctions*> > bfv(2);
  bfv[0]=bf; bfv[1]=bf2;
  std::vector<std::vector<Value*> > argsv(2);
  argsv[0]=args; argsv[1]=args;
  CoeffsVector mcoeffsV1 = CoeffsVector("mcoeffsV1",argsv,bfv,comm,true);
  mcoeffsV1.randomizeValuesGaussian(3);
  mcoeffsV1.writeToFile("mcoeffsV1.data",true);

  CoeffsMatrix mcoeffsM1 = CoeffsMatrix("mcoeffsM1",&mcoeffsV1,comm);
  mcoeffsM1.randomizeValuesGaussian(1);
  mcoeffsM1.writeToFile("mcoeffsM1.data");

  CoeffsVector coeffsV1 = CoeffsVector("coeffsV1",args,bf,comm,true);
  CoeffsVector coeffsV2 = CoeffsVector("coeffsV2",args,bf2,comm,true);


  std::vector<double> values(coeffsV1.numberOfCoeffs(),1.0);

  coeffsV1.randomizeValuesGaussian(1);
  coeffsV1.writeToFile("c1.1.data",true);
  coeffsV2 = coeffsV1-values;
  coeffsV2.writeToFile("c2.1.data",true);
  coeffsV2 = values-coeffsV1;
  coeffsV2.writeToFile("c2.2.data",true);
  coeffsV2 = coeffsV1-1.0;
  coeffsV2.writeToFile("c2.3.data",true);
  coeffsV2 = 1.0-coeffsV1;
  coeffsV2.writeToFile("c2.4.data",true);



  CoeffsMatrix coeffsM2 = CoeffsMatrix("Hessian1",&coeffsV1,comm,false);
  coeffsM2.randomizeValuesGaussian(1);
  coeffsM2.writeToFile("coeffsM2.data");

  std::vector<double> values2(coeffsM2.getSize(),1.0);

  CoeffsMatrix coeffsM3 = CoeffsMatrix("Hessian2",&coeffsV1,comm,false);
  coeffsM3.randomizeValuesGaussian(1);
  coeffsM3 = values2;
  coeffsM3.writeToFile("coeffsM3.data");

  CoeffsVector cV1 = CoeffsVector("cV1",args,bf,comm,true);
  CoeffsVector cV2 = CoeffsVector("cV2",args,bf,comm,true);
  cV1.randomizeValuesGaussian(1);
  cV2.randomizeValuesGaussian(2);
  cV1.writeToFile("cV1.1.data");
  cV2.writeToFile("cV2.1.data");
  std::vector<CoeffsVector*> cVv(2);
  cVv[0]=&cV1; cVv[1]=&cV2;
  CoeffsVector::averageVectors(cVv);
  cV1.writeToFile("cV1.2.data");
  cV1.writeToFile("cV2.2.data");




  CoeffsMatrix cM1 = CoeffsMatrix("cM1",args,bf,comm,true,true);
  CoeffsMatrix cM2 = CoeffsMatrix("cM2",args,bf,comm,true,true);
  cM1.randomizeValuesGaussian(1);
  cM2.randomizeValuesGaussian(2);
  cM1.writeToFile("cM1.1.data");
  cM2.writeToFile("cM2.1.data");
  std::vector<CoeffsMatrix*> cMv(2);
  cMv[0]=&cM1; cMv[1]=&cM2;
  CoeffsMatrix::averageMatrices(cMv);
  cM1.writeToFile("cM1.2.data");
  cM2.writeToFile("cM2.2.data");






  std::vector<std::string> min(2);
  std::vector<std::string> max(2);
  std::vector<unsigned int> nbins(2);
  min[0]="-4.0";
  min[1]="-4.0";
  max[0]="4.0";
  max[1]="4.0";
  nbins[0]=200;
  nbins[1]=200;

  std::string keywords = "GAUSSIAN CENTER0=-2.0,0.0 SIGMA0=0.5,0.5 CENTER1=+2.0,0.0 SIGMA1=0.5,0.5 WEIGHTS=1.0,10.0";
  TargetDistributionBase::writeDistributionToFile("dist",keywords,min,max,nbins);

  keywords = "LINEAR_COMBINATION DISTRIBUTION0={GAUSSIAN CENTER0=-2.0,0.0 SIGMA0=0.5,0.5} DISTRIBUTION1={GAUSSIAN CENTER0=+2.0,0.0 SIGMA0=0.5,0.5} WEIGHTS=1.0,10.0,2.0 DISTRIBUTION2={UNIFORM MINIMA=-2.0,-2.0 MAXIMA=2.0,1.0}";
  TargetDistributionBase::writeDistributionToFile("dist2",keywords,min,max,nbins);

  // keywords = "UNIFORM MINIMA=-2.0,-2.0 MAXIMA=2.0,1.0";
  //
  // keywords = "GRID FILE=dist ARGS=arg1,arg2 LABEL=GAUSSIAN";
  // TargetDistributionBase::writeDistributionToFile("dist3",keywords,min,max,nbins);
  //
  // keywords = "GRID FILE=dist ARGS=arg1,arg2 LABEL=GAUSSIAN NORMALIZE";
  // TargetDistributionBase::writeDistributionToFile("dist4",keywords,min,max,nbins);





}

void TestBFs::calculate(){

  std::vector<double> values(bf_pointer->getNumberOfBasisFunctions());
  std::vector<double> derivs(bf_pointer->getNumberOfBasisFunctions());
  bool inside_interval=true;
  double arg=getArgument(0);
  double argT=0.0;
  bf_pointer->getAllValues(arg,argT,inside_interval,values,derivs);

  getPntrToComponent("value")->set(values[bf_order_]);
  getPntrToComponent("deriv")->set(derivs[bf_order_]);
  getPntrToComponent("argT")->set(argT);
  getPntrToComponent("arg")->set(arg);
  if(inside_interval){getPntrToComponent("inside")->set(1.0);}
  else{getPntrToComponent("inside")->set(0.0);}
}

}
}
