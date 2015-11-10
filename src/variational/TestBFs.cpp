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
#include "../function/Function.h"
#include "core/ActionRegister.h"
#include "core/ActionSet.h"
#include "core/PlumedMain.h"
#include "BasisFunctions.h"
#include "CoeffsVector.h"
#include "CoeffsMatrix.h"
#include "tools/File.h"
#include "LinearBiasExpansion.h"
#include "tools/Communicator.h"
#include "TargetDistributionBase.h"


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
  CoeffsVector* coeffs2;
  LinearBiasExpansion* bias_expansion;
  unsigned int bf_order_;
public:
  TestBFs(const ActionOptions&);
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
  // bf_pointer->printInfo();
  parse("N",bf_order_);
  addComponent("value"); componentIsNotPeriodic("value");
  addComponent("deriv"); componentIsNotPeriodic("deriv");
  addComponent("arg"); componentIsNotPeriodic("arg");
  addComponent("argT"); componentIsNotPeriodic("argT");
  addComponent("inside"); componentIsNotPeriodic("inside");
  checkRead();
  log.printf("  using the %d order basis function from the %s basis set\n",bf_order_,basisset_label.c_str());

  std::vector<BasisFunctions*> bf; bf.resize(2); bf[0]=bf_pointer; bf[1]=bf_pointer;
  std::vector<Value*> args; args.resize(2); args[0]=getArguments()[0]; args[1]=getArguments()[1];
  CoeffsMatrix coeffsM = CoeffsMatrix("coeffsM",args,bf,comm);
  coeffsM.randomizeValuesGaussian(1);
  coeffsM.writeToFile("coeffsM.data");

  std::vector<unsigned int> ind1 = coeffsM.getIndices(4);
  std::vector<unsigned int> ind2 = coeffsM.getIndices(5);
  coeffsM(ind1,ind1)=4.0;
  coeffsM(ind1,ind1)+=23.0;

  coeffsM.writeToFile("coeffsM.1.data");
  coeffsM.sumMPI();
  coeffsM.writeToFile("coeffsM.2.data");

  // unsigned int rank comm.Rank();





  CoeffsVector coeffsV1 = CoeffsVector("coeffs",args,bf,comm,true);
  CoeffsVector coeffsV2 = CoeffsVector("coeffs",args,bf,comm,true);

  std::vector<double> vec( coeffsV1.getSize() );
  for (unsigned int i = 0; i < vec.size(); i++) {
    vec[i]=i*i;
  }
  coeffsV1 = 3.0;
  coeffsV2 = 100.0;
  CoeffsVector coeffsV3 = coeffsM*coeffsV1;
  CoeffsVector coeffsV4 = coeffsV1*coeffsM;

  // CoeffsVector coeffsV4 = 1.0-coeffsV1;
  CoeffsVector coeffsV5 = coeffsV1+1.0;
  CoeffsVector coeffsV6 = coeffsV1;
  CoeffsVector coeffsV7(coeffsV1);
  coeffsV7 = vec;
  coeffsV1.clear();

  for(unsigned int i=comm.Get_rank(); i<coeffsV1.getSize(); i+=comm.Get_size()){
    coeffsV1(i)=3.0;
  }
  coeffsV1.writeToFile("coeffsV1.1.data");
  coeffsV1.sumMPI();
  coeffsV1.writeToFile("coeffsV1.2.data");

  std::vector<std::string> min(2);
  std::vector<std::string> max(2);
  std::vector<unsigned int> nbins(2);

  min[0]="-4.0";
  min[1]="-4.0";
  max[0]="4.0";
  max[1]="4.0";
  nbins[0]=200;
  nbins[1]=200;

  std::string keywords = "GAUSSIAN CENTER0=-0.0,0.0 SIGMA0=0.5,0.5 CORRELATION0=-0.89";
  TargetDistributionBase::writeDistributionToFile("dist",keywords,min,max,nbins);

  keywords = "UNIFORM MINIMA=-2.0,-2.0 MAXIMA=2.0,1.0";
  TargetDistributionBase::writeDistributionToFile("dist2",keywords,min,max,nbins);

  keywords = "GRID FILE=dist ARGS=arg1,arg2 LABEL=GAUSSIAN";
  TargetDistributionBase::writeDistributionToFile("dist3",keywords,min,max,nbins);

  keywords = "GRID FILE=dist ARGS=arg1,arg2 LABEL=GAUSSIAN NORMALIZE";
  TargetDistributionBase::writeDistributionToFile("dist4",keywords,min,max,nbins);



  // std::vector<CoeffsVector> d1;
  // d1.push_back(coeffsV1);
  // CoeffsVector::writeToFile("test.data",d1,comm,true);

  // coeffsV2.writeToFile("coeffsV2.data");
  // coeffsV3.writeToFile("coeffsV3.data");
  // coeffsV4.writeToFile("coeffsV4.data");
  // coeffsV5.writeToFile("coeffsV5.data");
  // coeffsV6.writeToFile("coeffsV6.data");
  // coeffsV7.writeToFile("coeffsV7.data");


    /*
  std::vector<BasisFunctions*> bf2; bf2.resize(1); bf2[0]=bf_pointer;
  std::vector<Value*> args2; args2.resize(1); args2[0]=getArguments()[0];
  CoeffsVector* coeffsV2 = new CoeffsVector("coeffs2",args2,bf2,true);
  // d1.push_back(*coeffsV2);
  (*coeffsV2).writeToFile("dd.data");

  CoeffsVector coeffsV_copy1 = CoeffsVector(*coeffsV);
  coeffsV_copy1.setValues(10.0);
  coeffsV_copy1.setDataLabel("aux_coeffs");
  d1.push_back(coeffsV_copy1);
  coeffsV_copy1.clear();
  coeffsV_copy1.writeToFile("coeffsV_copy1.before.data");




  CoeffsVector coeffsV_copy2 = 2*CoeffsVector(coeffsV_copy1);
  coeffsV_copy2.setDataLabel("aux_coeffs2");
  coeffsV_copy2.setValues(2.0);
  d1.push_back(coeffsV_copy2);
  coeffsV_copy2.clear();
  CoeffsV_copy2.writeToFile("coeffsV_copy2.before.data");



  d1.clear();
  d1.push_back(*coeffsV);
  CoeffsVector::writeToFile("test.data",d1,true,true);

  coeffsV->setValues(3e56);
  coeffsV_copy1.setValues(3.0);
  coeffsV_copy2.randomizeValuesGaussian(414);

  d1.clear();
  d1.push_back(*coeffsV);
  d1.push_back(coeffsV_copy1);
  d1.push_back(coeffsV_copy2);
  CoeffsVector::writeToFile("test.data",d1,true,true);
  */


  /*
  std::vector<std::string> bf1;
  bf1.push_back("BF_FOURIER");
  bf1.push_back("ORDER=10");
  bf1.push_back("LABEL=bf2");
  bf1.push_back("INTERVAL_MIN=-pi");
  bf1.push_back("INTERVAL_MAX=pi");
  plumed.readInputWords(bf1);
  // BasisFunctions* bf_pointer2=plumed.getActionSet().selectWithLabel<BasisFunctions*>("bf2");

  std::vector<BasisFunctions*> bf; bf.resize(2); bf[0]=bf_pointer; bf[1]=bf_pointer;
  std::vector<Value*> args; args.resize(2); args[0]=getArguments()[0]; args[1]=getArguments()[1];
  bias_expansion = new LinearBiasExpansion("bla",args,bf,comm);
  std::vector<unsigned int> nbins(2,300);
  bias_expansion->setupGrid(nbins);

  coeffs2 = new CoeffsVector("Test",args,bf,true,true);
  // for(unsigned int i=0;i<coeffs2->getSize();i++){coeffs2->setValue(i,1.0*i*i*i);}
  coeffs2->randomizeValuesGaussian();
  coeffs2->setCounter(100);
  coeffs2->writeToFile("TEST2.data",true,true);
  log.printf("Min:  %f\n",coeffs2->getMinValue());
  log.printf("Max:  %f\n",coeffs2->getMaxValue());
  log.printf("Norm: %f\n",coeffs2->getNorm());
  coeffs2->normalizeCoeffs();
  coeffs2->writeToFile("TEST3.data",true,true);
  coeffs2->writeToFile("TEST4.data",true,true);
  coeffs2->writeToFile("TEST5.data",true,true);


  bias_expansion->updateBiasGrid();
  bias_expansion->writeBiasGridToFile("bias.data",false);
  bias_expansion->writeBiasGridToFile("bias2.data",false);
  bias_expansion->writeBiasGridToFile("bias2.data",true);

  TargetDistribution1DimBase::writeDistributionToFile("dist","GAUSSIAN CENTER=-2.0,2.0 SIGMA=0.5,0.5 WEIGHT=1,4 DO_NOT_NORMALIZE",-4.0,4.0,200);


  */

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
