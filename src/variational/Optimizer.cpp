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
#include "Optimizer.h"
#include "Coeffs.h"
#include "VariationalBias.h"


namespace PLMD{

Optimizer::Optimizer(const ActionOptions&ao):
Action(ao),
ActionPilot(ao),
ActionWithValue(ao),
has_been_set(false),
usehessian_(false),
description_("Undefined"),
type_("Undefined")
{
 parse("STEP_SIZE",step_size_);
}

void Optimizer::registerKeywords( Keywords& keys ){
  Action::registerKeywords(keys);
  ActionPilot::registerKeywords(keys);
  ActionWithValue::registerKeywords(keys);
  keys.add("compulsory","STEP_SIZE","the step size used for the optimization");
}

void Optimizer::apply(){}

void Optimizer::calculate(){}

void Optimizer::linkCoeffs(Coeffs* bias_coeffs_in, Coeffs* gradient_in, Coeffs* hessian_diag_in)
{
 if(usehessian_==false){plumed_merror("link to hessian is not needed");}
 bias_coeffs = bias_coeffs_in;
 gradient = gradient_in;
 hessian_diag = hessian_diag_in;
}

void Optimizer::linkCoeffs(Coeffs* bias_coeffs_in, Coeffs* gradient_in)
{
 if(usehessian_==true){plumed_merror("link to hessian also needed");}
 bias_coeffs = bias_coeffs_in;
 gradient = gradient_in;
}

void Optimizer::linkBias(VariationalBias* bias_ptr_in){bias_ptr=bias_ptr_in;}


}


