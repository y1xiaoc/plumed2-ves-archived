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
#include "core/ActionRegister.h"
#include "BasisFunctions.h"

namespace PLMD{

class BF_Cosine : public BasisFunctions {
  virtual void setupLabels();
  virtual void setupUniformIntegrals();
public:
  static void registerKeywords(Keywords&);
  explicit BF_Cosine(const ActionOptions&);
  double getValue(const double, const unsigned int, double&, bool&) const;
  void getAllValues(const double, double&, bool&, std::vector<double>&, std::vector<double>&) const;
};


PLUMED_REGISTER_ACTION(BF_Cosine,"BF_COSINE")


void BF_Cosine::registerKeywords(Keywords& keys){
  BasisFunctions::registerKeywords(keys);
}


BF_Cosine::BF_Cosine(const ActionOptions&ao):
PLUMED_BASISFUNCTIONS_INIT(ao)
{
  setNumberOfBasisFunctions(getOrder()+1);
  setIntrinsicInterval("-pi","+pi");
  setPeriodic();
  setIntervalBounded();
  setType("trigonometric_cos");
  setDescription("Cosine");
  setupBF();
  checkRead();
}


double BF_Cosine::getValue(const double arg, const unsigned int n, double& argT, bool& inside_range) const {
  plumed_massert(n<numberOfBasisFunctions(),"getValue: n is outside range of the defined order of the basis set");
  inside_range=true;
  argT=translateArgument(arg, inside_range);
  double value=0.0;
  if(n == 0){
    value=1.0;
  }
  else {
    double k = n;
    value=cos(k*argT);
  }
  return value;
}


void BF_Cosine::getAllValues(const double arg, double& argT, bool& inside_range, std::vector<double>& values, std::vector<double>& derivs) const {
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
    values[i] = cos_tmp;
    derivs[i] = -io*sin_tmp*intervalDerivf();
  }
  if(!inside_range){
    for(unsigned int i=0;i<derivs.size();i++){derivs[i]=0.0;}
  }
}


void BF_Cosine::setupLabels() {
  setLabel(0,"1");
  for(unsigned int i=1; i < getOrder()+1;i++){
    std::string is; Tools::convert(i,is);
    setLabel(i,"cos("+is+"*s)");
  }
}


void BF_Cosine::setupUniformIntegrals() {
  setAllUniformIntegralsToZero();
  setUniformIntegral(0,1.0);
}


}
