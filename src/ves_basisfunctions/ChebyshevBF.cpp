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

class ChebyshevBF : public BasisFunctions {
  virtual void setupUniformIntegrals();
public:
  static void registerKeywords(Keywords&);
  explicit ChebyshevBF(const ActionOptions&);
  double getValue(const double, const unsigned int, double&, bool&) const;
  void getAllValues(const double, double&, bool&, std::vector<double>&, std::vector<double>&) const;
};


PLUMED_REGISTER_ACTION(ChebyshevBF,"BF_CHEBYSHEV")


void ChebyshevBF::registerKeywords(Keywords& keys){
 BasisFunctions::registerKeywords(keys);
}

ChebyshevBF::ChebyshevBF(const ActionOptions&ao):
 PLUMED_BASISFUNCTIONS_INIT(ao)
{
  setNumberOfBasisFunctions(getOrder()+1);
  setIntrinsicInterval(-1.0,+1.0);
  setNonPeriodic();
  setIntervalBounded();
  setType("chebyshev-1st-kind");
  setDescription("Chebyshev polynomials of the first kind");
  setLabelPrefix("T");
  setupBF();
}


double ChebyshevBF::getValue(const double arg, const unsigned int n, double& argT, bool& inside_range) const {
  plumed_massert(n<numberOfBasisFunctions(),"getValue: n is outside range of the defined order of the basis set");
  inside_range=true;
  std::vector<double> tmp_values(numberOfBasisFunctions());
  std::vector<double> tmp_derivs(numberOfBasisFunctions());
  getAllValues(arg, argT, inside_range, tmp_values, tmp_derivs);
  return tmp_values[n];
}


void ChebyshevBF::getAllValues(const double arg, double& argT, bool& inside_range, std::vector<double>& values, std::vector<double>& derivs) const {
  // plumed_assert(values.size()==numberOfBasisFunctions());
  // plumed_assert(derivs.size()==numberOfBasisFunctions());
  inside_range=true;
  argT=translateArgument(arg, inside_range);
  std::vector<double> derivsT(derivs.size());
  //
  values[0]=1.0;
  derivsT[0]=0.0;
  derivs[0]=0.0;
  values[1]=argT;
  derivsT[1]=1.0;
  derivs[1]=intervalDerivf();
  for(unsigned int i=1; i < getOrder();i++){
    values[i+1]  = 2.0*argT*values[i]-values[i-1];
    derivsT[i+1] = 2.0*values[i]+2.0*argT*derivsT[i]-derivsT[i-1];
    derivs[i+1]  = intervalDerivf()*derivsT[i+1];
  }
  if(!inside_range){for(unsigned int i=0;i<derivs.size();i++){derivs[i]=0.0;}}
}


void ChebyshevBF::setupUniformIntegrals() {
  for(unsigned int i=0; i<numberOfBasisFunctions(); i++){
    double io = i;
    double value = 0.0;
    if(i % 2 == 0){
      value = -2.0/( pow(io,2.0)-1.0)*0.5;
    }
    setUniformIntegral(i,value);
  }
}


}
