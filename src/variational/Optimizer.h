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
#ifndef __PLUMED_variational_Optimizer_h
#define __PLUMED_variational_Optimizer_h

#include <vector>
#include <string>
#include <cmath>
#include "core/ActionWithValue.h"

#define PLUMED_BASISFUNCTIONS_INIT(ao) Action(ao),Optimizer(ao)

namespace PLMD{

/**
\ingroup INHERIT
Abstract base class for implenting new optimization methods
*/

class Coeffs;

class Optimizer :
 public ActionWithValue
{
private:

protected:
  bool has_been_set;
  bool usehessian_;
  // description of optimizer
  std::string description_;
  // the type of optimizer
  std::string type_;
  double step_size_;
  Coeffs* bias_coeffs;
  Coeffs* gradient;
  Coeffs* hessian_diag;
public:
  static void registerKeywords(Keywords&);
  Optimizer(const ActionOptions&ao);
  bool hasBeenSet();
  std::string getType();
  std::string getDescription();
  //
  void apply();
  void calculate();
  // 
  void linkCoeffs(Coeffs*, Coeffs*, Coeffs*);
  void linkCoeffs(Coeffs*, Coeffs*);

};

inline
bool Optimizer::hasBeenSet(){return has_been_set;}

inline
std::string Optimizer::getType(){return type_;}

inline
std::string Optimizer::getDescription(){return description_;}

}

#endif

