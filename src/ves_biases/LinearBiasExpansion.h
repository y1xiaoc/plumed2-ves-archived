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
#ifndef __PLUMED_ves_biases_LinearBiasExpansion_h
#define __PLUMED_ves_biases_LinearBiasExpansion_h

#include <vector>
#include <string>

namespace PLMD {

class Keywords;
class CoeffsVector;
class BasisFunctions;
class Value;
class Communicator;
class Grid;

class LinearBiasExpansion{
  Communicator& mycomm;
  bool serial_;
  std::string bias_label_;
  CoeffsVector* bias_coeffs;
  CoeffsVector* wt_coeffs;
  CoeffsVector* basisf_norm;
  Grid* bias_grid;
  Grid* fes_grid;
  Grid* ps_grid;
  std::vector<Value*> args_;
  std::vector<BasisFunctions*> basisf_;
  unsigned int ncv_;
  std::vector<unsigned int> num_bf_;
  //
 public:
  static void registerKeywords( Keywords& keys );
  // Constructor
  explicit LinearBiasExpansion(const std::string,
                      std::vector<Value*>,
                      std::vector<BasisFunctions*>,
                      Communicator &cc);
  //
  std::vector<Value*> getPointerToArguments() const ;
  std::vector<BasisFunctions*> getPointerToBasisFunctions() const ;
  CoeffsVector* getPointerToBiasCoeffs() const ;
  Grid* getPointerToBiasGrid() const ;
  unsigned int getNumberOfArguments() const ;
  std::vector<unsigned int> getNumberOfBasisFunctions() const ;
  unsigned int getNumberOfCoeffs() const ;
  // Grid stuff
  void setupGrid(const std::vector<unsigned int>&);
  void updateBiasGrid();
  void writeBiasGridToFile(const std::string, const bool);
  // calculate bias and derivatives
  double getBiasAndDerivatives(const std::vector<double>&, std::vector<double>& derivatives);
};

}

#endif
