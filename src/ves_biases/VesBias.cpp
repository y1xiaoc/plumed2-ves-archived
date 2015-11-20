/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2011-2015 The plumed team
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
#include "VesBias.h"
#include "ves_basisfunctions/BasisFunctions.h"
#include "ves_tools/CoeffsVector.h"
#include "ves_tools/CoeffsMatrix.h"


namespace PLMD{
namespace bias{

VesBias::VesBias(const ActionOptions&ao): Bias(ao){}

VesBias::~VesBias(){
  delete coeffs_ptr;
  delete gradient_ptr;
  delete hessian_ptr;
}


void VesBias::registerKeywords( Keywords& keys ){
  Bias::registerKeywords(keys);
}


void VesBias::initializeCoeffs(const std::vector<std::string>& dimension_labels,const std::vector<unsigned int>& indices_shape) {
  coeffs_ptr = new CoeffsVector("coeffs",dimension_labels,indices_shape,comm,true);
  initializeGradientAndHessian();
}

void VesBias::initializeCoeffs(std::vector<Value*>& args,std::vector<BasisFunctions*>& basisf) {
  coeffs_ptr = new CoeffsVector("coeffs",args,basisf,comm,true);
  initializeGradientAndHessian();
}


void VesBias::initializeGradientAndHessian() {
  //
  gradient_ptr = new CoeffsVector(*coeffs_ptr);
  gradient_ptr->setLabel("gradient");
  gradient_ptr->setDataLabel("gradient");
  hessian_ptr = new CoeffsMatrix("hessian",coeffs_ptr,comm);

}


void VesBias::updateGradientAndHessian(){}


void VesBias::clearGradientAndHessian(){}

}
}
