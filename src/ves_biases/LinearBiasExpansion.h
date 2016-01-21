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

class Action;
class Keywords;
class Value;
class Communicator;
class Grid;
class CoeffsVector;
class BasisFunctions;
namespace bias{
  class VesBias;


class LinearBiasExpansion{
private:
  std::string label_;
  //
  Action* action_pntr;
  bias::VesBias* vesbias_pntr;
  Communicator& mycomm;
  bool serial_;
  //
  std::vector<Value*> args_pntrs;
  unsigned int nargs_;
  //
  std::vector<BasisFunctions*> basisf_pntrs;
  std::vector<unsigned int> nbasisf_;
  //
  CoeffsVector* bias_coeffs_pntr;
  size_t ncoeffs_;
  CoeffsVector* coeffderivs_aver_ps_pntr;
  CoeffsVector* fes_wt_coeffs_pntr;
  //
  double biasf_;
  double invbiasf_;
  //
  Grid* bias_grid_pntr;
  Grid* fes_grid_pntr;
  Grid* ps_grid_pntr;
  //
 public:
  static void registerKeywords( Keywords& keys );
  // Constructor
  explicit LinearBiasExpansion(
    const std::string&,
    Communicator &cc,
    std::vector<Value*>,
    std::vector<BasisFunctions*>,
    CoeffsVector* bias_coeffs_pntr_in=NULL);
  //
  ~LinearBiasExpansion();
  //
  std::vector<Value*> getPntrsToArguments() const ;
  std::vector<BasisFunctions*> getPntrsToBasisFunctions() const ;
  CoeffsVector* getPntrToBiasCoeffs() const ;
  Grid* getPntrToBiasGrid() const ;
  //
  unsigned int getNumberOfArguments() const ;
  std::vector<unsigned int> getNumberOfBasisFunctions() const ;
  size_t getNumberOfCoeffs() const ;
  //
  CoeffsVector& BiasCoeffs() const;
  CoeffsVector& FesWTCoeffs() const;
  //
  void setSerial() {serial_=true;}
  void setParallel() {serial_=false;}
  //
  void linkVesBias(bias::VesBias*);
  void linkAction(Action*);
  // calculate bias and derivatives
  double getBiasAndForces(const std::vector<double>&, std::vector<double>&);
  double getBias(const std::vector<double>&);
  // Grid stuff
  void setupGrid(const std::vector<unsigned int>&, const bool usederiv=false);
  void updateBiasGrid();
  void writeBiasGridToFile(const std::string&, const bool);
  // Well-Tempered p(s) stuff
  void setupWellTempered(const double, const std::vector<unsigned int>&);
  void updateWellTemperedFESCoeffs();
};

inline
std::vector<Value*> LinearBiasExpansion::getPntrsToArguments() const {return args_pntrs;}

inline
std::vector<BasisFunctions*> LinearBiasExpansion::getPntrsToBasisFunctions() const {return basisf_pntrs;}

inline
CoeffsVector* LinearBiasExpansion::getPntrToBiasCoeffs() const {return bias_coeffs_pntr;}

inline
Grid* LinearBiasExpansion::getPntrToBiasGrid() const {return bias_grid_pntr;}

inline
unsigned int LinearBiasExpansion::getNumberOfArguments() const {return nargs_;}

inline
std::vector<unsigned int> LinearBiasExpansion::getNumberOfBasisFunctions() const {return nbasisf_;}

inline
size_t LinearBiasExpansion::getNumberOfCoeffs() const {return ncoeffs_;}

inline
CoeffsVector& LinearBiasExpansion::BiasCoeffs() const {return *bias_coeffs_pntr;}

inline
CoeffsVector& LinearBiasExpansion::FesWTCoeffs() const {return *fes_wt_coeffs_pntr;}

}

}

#endif
