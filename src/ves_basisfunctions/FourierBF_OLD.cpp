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
#include "BasisFunctions.h"

#include "core/ActionRegister.h"

namespace PLMD{

class FourierBF_OLD : public BasisFunctions {
  virtual void setupLabels();
  virtual void setupUniformIntegrals();
public:
  static void registerKeywords(Keywords&);
  explicit FourierBF_OLD(const ActionOptions&);
  double getValue(const double, const unsigned int, double&, bool&) const;
  void getAllValues(const double, double&, bool&, std::vector<double>&, std::vector<double>&) const;
};


PLUMED_REGISTER_ACTION(FourierBF_OLD,"BF_FOURIER_OLD")


void FourierBF_OLD::registerKeywords(Keywords& keys){
  BasisFunctions::registerKeywords(keys);
}


FourierBF_OLD::FourierBF_OLD(const ActionOptions&ao):
PLUMED_BASISFUNCTIONS_INIT(ao)
{
  setNumberOfBasisFunctions(2*getOrder()+1);
  setIntrinsicInterval(-pi,+pi);
  setPeriodic();
  setIntervalBounded();
  setType("trigonometric_cos-sin");
  setDescription("Trigonometric (cos/sin)");
  setupBF();
}


double FourierBF_OLD::getValue(const double arg, const unsigned int n, double& argT, bool& inside_range) const {
  plumed_massert(n<numberOfBasisFunctions(),"getValue: n is outside range of the defined order of the basis set");
  inside_range=true;
  argT=translateArgument(arg, inside_range);
  double value=0.0;
  double n_cos = n;
  double n_sin = n - getOrder();
  if(n==0){value=1.0;}
  else if(n <= getOrder()){value=cos(n_cos*argT);}
  else if(n <= 2*getOrder()){value=sin(n_sin*argT);}
  return value;
}


void FourierBF_OLD::getAllValues(const double arg, double& argT, bool& inside_range, std::vector<double>& values, std::vector<double>& derivs) const {
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
    values[i+getOrder()] = sin_tmp;
    derivs[i+getOrder()] = io*cos_tmp*intervalDerivf();
  }
  if(!inside_range){for(unsigned int i=0;i<derivs.size();i++){derivs[i]=0.0;}}
}


void FourierBF_OLD::setupLabels() {
  setLabel(0,"1");
  for(unsigned int i=1; i < getOrder()+1;i++){
    std::string is; Tools::convert(i,is);
    setLabel(i,"cos("+is+"*s)");
    setLabel(i+getOrder(),"sin("+is+"*s)");
  }
}


void FourierBF_OLD::setupUniformIntegrals() {
  setAllUniformIntegralsToZero();
  setUniformIntegral(0,1.0);
}


}
