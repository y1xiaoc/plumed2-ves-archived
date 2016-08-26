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

class LegendreBF : public BasisFunctions {
  bool scaled_;
  virtual void setupUniformIntegrals();
public:
  static void registerKeywords(Keywords&);
  explicit LegendreBF(const ActionOptions&);
  double getValue(const double, const unsigned int, double&, bool&) const;
  void getAllValues(const double, double&, bool&, std::vector<double>&, std::vector<double>&) const;
};


PLUMED_REGISTER_ACTION(LegendreBF,"BF_LEGENDRE")


void LegendreBF::registerKeywords(Keywords& keys){
 BasisFunctions::registerKeywords(keys);
 keys.addFlag("SCALED",false,"scale the polynomials such that they are orthonormal to 1");
}

LegendreBF::LegendreBF(const ActionOptions&ao):
 PLUMED_BASISFUNCTIONS_INIT(ao),
 scaled_(false)
{
  parseFlag("SCALED",scaled_); addKeywordToList("SCALED",scaled_);
  setNumberOfBasisFunctions(getOrder()+1);
  setIntrinsicInterval("-1.0","+1.0");
  setNonPeriodic();
  setIntervalBounded();
  setType("Legendre");
  setDescription("Legendre polynomials");
  setLabelPrefix("L");
  setupBF();
  checkRead();
}


double LegendreBF::getValue(const double arg, const unsigned int n, double& argT, bool& inside_range) const {
  plumed_massert(n<numberOfBasisFunctions(),"getValue: n is outside range of the defined order of the basis set");
  inside_range=true;
  std::vector<double> tmp_values(numberOfBasisFunctions());
  std::vector<double> tmp_derivs(numberOfBasisFunctions());
  getAllValues(arg, argT, inside_range, tmp_values, tmp_derivs);
  return tmp_values[n];
}


void LegendreBF::getAllValues(const double arg, double& argT, bool& inside_range, std::vector<double>& values, std::vector<double>& derivs) const {
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
    double io = static_cast<double>(i);
    values[i+1]  = ((2.0*io+1.0)/(io+1.0))*argT*values[i] - (io/(io+1.0))*values[i-1];
    derivsT[i+1] = ((2.0*io+1.0)/(io+1.0))*(values[i]+argT*derivsT[i])-(io/(io+1.0))*derivsT[i-1];
    derivs[i+1]  = intervalDerivf()*derivsT[i+1];
  }
  if(scaled_){
    for(unsigned int i=0; i<values.size(); i++){
      double io = static_cast<double>(i);
      double sf = sqrt(io+0.5);
      values[i] *= sf;
      derivs[i] *= sf;
    }
  }
  if(!inside_range){for(unsigned int i=0;i<derivs.size();i++){derivs[i]=0.0;}}
}


void LegendreBF::setupUniformIntegrals() {
  setAllUniformIntegralsToZero();
  setUniformIntegral(0,1.0);
}


}
