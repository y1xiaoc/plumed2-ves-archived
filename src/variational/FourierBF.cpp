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
namespace BasisFunctions{

class FourierBF : public BasisFunctions{
 virtual void setupDescription();
 virtual void setupBFIntegrals();
public:
 FourierBF(const ActionOptions&);
 double getValue(const double, const unsigned int, double&, bool&);
 void getAllValues(const double, double&, bool&, std::vector<double>&, std::vector<double>&);
};

PLUMED_REGISTER_ACTION(FourierBF,"BF_FOURIER")

FourierBF::FourierBF(const ActionOptions&ao):
PLUMED_BASISFUNCTIONS_INIT(ao)
{
 nbasis_ = 2*norder_+1;
 interval_default_min_=-pi;
 interval_default_max_=+pi;
 periodic_=true;
 interval_bounded_=true;
 type_="Trigonometric (cos/sin)";
 setupBF();
 printInfo();
}

double FourierBF::getValue(const double arg, const unsigned int n, double& argT, bool& inside_range)
{
 if(n>=nbasis_){error("getValue: n is outside range of the defined order of the basis set");}
 inside_range=true;
 argT=translateArgument(arg, inside_range);
 double value=0.0;
 double n_cos = n;
 double n_sin = n-norder_;
 if(n==0){value=1.0;}
 else if(n <= norder_){value=cos(n_cos*argT);}
 else if(n <= 2*norder_){value=sin(n_sin*argT);}
 return value;
}

void FourierBF::getAllValues(const double arg, double& argT, bool& inside_range, std::vector<double>& values, std::vector<double>& derivs)
{
 if(values.size()!=nbasis_ || derivs.size()!=nbasis_){error("getAllValues: wrong size of values or derivs vectors");}
 inside_range=true;
 argT=translateArgument(arg, inside_range);
 double cos_tmp=0; 
 double sin_tmp=0; 
 values[0]=1.0;
 derivs[0]=0.0;
 for(unsigned int i=1; i < norder_+1;i++)
 {
  double io = i;
  cos_tmp = cos(io*argT);
  sin_tmp = sin(io*argT);
  values[i] = cos_tmp;
  derivs[i] = -io*sin_tmp*argT_derivf_;
  values[i+norder_] = sin_tmp;
  derivs[i+norder_] = io*cos_tmp*argT_derivf_;
 }
 if(!inside_range){for(unsigned int i=0;i<derivs.size();i++){derivs[i]=0.0;}} 
}

void FourierBF::setupDescription()
{
 bf_description_.resize(nbasis_);
 bf_description_[0]="1";
 for(unsigned int i=1; i < norder_+1;i++)
 {
  std::string is; Tools::convert(i,is);
  bf_description_[i]="cos("+is+"*s)";
  bf_description_[i+norder_]="sin("+is+"*s)";
 }
}

void FourierBF::setupBFIntegrals()
{
 bf_integrals_.assign(nbasis_,0.0);
 bf_integrals_[0]=1.0;
}


}
}


