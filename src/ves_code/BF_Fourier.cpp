/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2015-2016 The ves-code team
   (see the PEOPLE-VES file at the root of the distribution for a list of names)

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

#include "BasisFunctions.h"

#include "core/ActionRegister.h"


namespace PLMD{

class BF_Fourier : public BasisFunctions {
  virtual void setupLabels();
  virtual void setupUniformIntegrals();
public:
  static void registerKeywords(Keywords&);
  explicit BF_Fourier(const ActionOptions&);
  double getValue(const double, const unsigned int, double&, bool&) const;
  void getAllValues(const double, double&, bool&, std::vector<double>&, std::vector<double>&) const;
};


PLUMED_REGISTER_ACTION(BF_Fourier,"BF_FOURIER")


void BF_Fourier::registerKeywords(Keywords& keys){
  BasisFunctions::registerKeywords(keys);
}


BF_Fourier::BF_Fourier(const ActionOptions&ao):
PLUMED_BASISFUNCTIONS_INIT(ao)
{
  setNumberOfBasisFunctions(2*getOrder()+1);
  setIntrinsicInterval("-pi","+pi");
  setPeriodic();
  setIntervalBounded();
  setType("trigonometric_cos-sin");
  setDescription("Trigonometric (cos/sin)");
  setupBF();
  checkRead();
}


double BF_Fourier::getValue(const double arg, const unsigned int n, double& argT, bool& inside_range) const {
  plumed_massert(n<numberOfBasisFunctions(),"getValue: n is outside range of the defined order of the basis set");
  inside_range=true;
  argT=translateArgument(arg, inside_range);
  double value=0.0;
  if(n == 0){
    value=1.0;
  }
  else if(n%2 == 1){
    double k = (n+1.0)/2.0;
    value=cos(k*argT);
  }
  else if(n%2 == 0){
    double k = n/2.0;
    value=sin(k*argT);
  }
  return value;
}


void BF_Fourier::getAllValues(const double arg, double& argT, bool& inside_range, std::vector<double>& values, std::vector<double>& derivs) const {
  // plumed_assert(values.size()==numberOfBasisFunctions());
  // plumed_assert(derivs.size()==numberOfBasisFunctions());
  inside_range=true;
  argT=translateArgument(arg, inside_range);
  values[0]=1.0;
  derivs[0]=0.0;
  for(unsigned int i=1; i < getOrder()+1;i++){
    double io = i;
    double cos_tmp = cos(io*argT);
    double sin_tmp = sin(io*argT);
    values[2*i-1] = cos_tmp;
    derivs[2*i-1] = -io*sin_tmp*intervalDerivf();
    values[2*i] = sin_tmp;
    derivs[2*i] = io*cos_tmp*intervalDerivf();
  }
  if(!inside_range){
    for(unsigned int i=0;i<derivs.size();i++){derivs[i]=0.0;}
  }
}


void BF_Fourier::setupLabels() {
  setLabel(0,"1");
  for(unsigned int i=1; i < getOrder()+1;i++){
    std::string is; Tools::convert(i,is);
    setLabel(2*i-1,"cos("+is+"*s)");
    setLabel(2*i,"sin("+is+"*s)");
  }
}


void BF_Fourier::setupUniformIntegrals() {
  setAllUniformIntegralsToZero();
  setUniformIntegral(0,1.0);
}


}
