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

class CubicBsplineBF : public BasisFunctions {
  double spacing_;
  double inv_spacing_;
  double spline(const double, double&) const;
public:
  static void registerKeywords( Keywords&);
  explicit CubicBsplineBF(const ActionOptions&);
  double getValue(const double, const unsigned int, double&, bool&) const;
  void getAllValues(const double, double&, bool&, std::vector<double>&, std::vector<double>&) const;
};


PLUMED_REGISTER_ACTION(CubicBsplineBF,"BF_CUBIC_B_SPLINES")


void CubicBsplineBF::registerKeywords(Keywords& keys){
  BasisFunctions::registerKeywords(keys);
}

CubicBsplineBF::CubicBsplineBF(const ActionOptions&ao):
PLUMED_BASISFUNCTIONS_INIT(ao)
{
  // plumed_merror("BF_CUBIC_B_SPLINES are not yet implented");
  setNumberOfBasisFunctions((getOrder()+3)+1);
  setIntrinsicInterval(intervalMin(),intervalMax());
  spacing_=(intervalMax()-intervalMin())/static_cast<double>(getOrder());
  inv_spacing_ = 1.0/spacing_;
  setNonPeriodic();
  setIntervalBounded();
  setType("splines_2nd-order");
  setDescription("Cubic B-splines (2nd order splines)");
  setLabelPrefix("S");
  setupBF();
}


double CubicBsplineBF::getValue(const double arg, const unsigned int n, double& argT, bool& inside_range) const {
  plumed_massert(n<numberOfBasisFunctions(),"getValue: n is outside range of the defined order of the basis set");
  inside_range=true;
  argT=translateArgument(arg, inside_range);
  //
  if(n==0){
    return 1.0;
  }
  else{
    double argx = ((argT-intervalMin())/spacing_) - (static_cast<double>(n)-2.0);
    double tmp_dbl=0.0;
    return spline(argx, tmp_dbl);
  }
}


void CubicBsplineBF::getAllValues(const double arg, double& argT, bool& inside_range, std::vector<double>& values, std::vector<double>& derivs) const {
  // plumed_assert(values.size()==numberOfBasisFunctions());
  // plumed_assert(derivs.size()==numberOfBasisFunctions());
  inside_range=true;
  argT=translateArgument(arg, inside_range);
  //
  values[0]=1.0;
  derivs[0]=0.0;
  //
  for(unsigned int i=1; i < getNumberOfBasisFunctions(); i++){
    double argx = ((argT-intervalMin())/spacing_) - (static_cast<double>(i)-2.0);
    values[i]  = spline(argx, derivs[i]);
    derivs[i]*=inv_spacing_;
  }
  if(!inside_range){for(unsigned int i=0;i<derivs.size();i++){derivs[i]=0.0;}}
}


double CubicBsplineBF::spline(const double arg, double& deriv) const {
  double value=0.0;
  double x=arg;
  // derivative of abs(x);
  double dx = 1.0;
  if(x < 0){
    x=-x;
    dx = -1.0;
  }
  //
  if(x > 2){
    value=0.0;
    deriv=0.0;
  }
  else if(x >= 1){
    value = ((2.0-x)*(2.0-x)*(2.0-x));
    deriv = dx*(-3.0*(2.0-x)*(2.0-x));
    // value=((2.0-x)*(2.0-x)*(2.0-x))/6.0;
    // deriv=-x*x*(2.0-x)*(2.0-x);
  }
  else{
    value = 4.0-6.0*x*x+3.0*x*x*x;
    deriv = dx*(-12.0*x+9.0*x*x);
    // value=x*x*x*0.5-x*x+2.0/3.0;
    // deriv=(3.0/2.0)*x*x-2.0*x;
  }
  deriv/=6.0;
  value/=6.0;
  return value;
}


}
