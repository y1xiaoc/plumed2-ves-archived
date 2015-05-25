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
#ifndef __PLUMED_variational_LinearBiasExpansion_h
#define __PLUMED_variational_LinearBiasExpansion_h

#include <vector>
#include <string>
// #include "Coeffs.h"

namespace PLMD {

class Keywords;
class Coeffs;
class BasisFunctions;
class Value;
class Communicator;

class LinearBiasExpansion{
  Communicator & comm;
  std::string bias_label_;
  bool serial_;
  Coeffs* bias_coeffs;
  std::vector<Value*> args_;
  std::vector<BasisFunctions*> basisf_;
  unsigned int ncv_;
  std::vector<unsigned int> num_bf_;
///
 public:
  static void registerKeywords( Keywords& keys );
/// Constructor
  LinearBiasExpansion(const std::string, 
                      std::vector<Value*>,
                      std::vector<BasisFunctions*>,
                      Communicator& cc);
  Coeffs* getPointerToCoeffs() const ;
  double getBiasAndDerivatives(const std::vector<double>&, double* derivatives=NULL);
};

}

#endif

