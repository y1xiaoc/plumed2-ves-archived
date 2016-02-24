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
#ifndef __PLUMED_ves_biases_LinearBasisSetExpansion_h
#define __PLUMED_ves_biases_LinearBasisSetExpansion_h

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
class TargetDistribution;

namespace bias{
  class VesBias;


class LinearBasisSetExpansion{
private:
  std::string label_;
  //
  Action* action_pntr_;
  bias::VesBias* vesbias_pntr_;
  Communicator& mycomm_;
  bool serial_;
  //
  std::vector<Value*> args_pntrs_;
  unsigned int nargs_;
  //
  std::vector<BasisFunctions*> basisf_pntrs_;
  std::vector<unsigned int> nbasisf_;
  //
  CoeffsVector* bias_coeffs_pntr_;
  size_t ncoeffs_;
  CoeffsVector* coeffderivs_aver_ps_pntr_;
  CoeffsVector* fes_wt_coeffs_pntr_;
  //
  double welltemp_biasf_;
  double inv_welltemp_biasf_;
  //
  //
  std::vector<unsigned int> grid_bins_;
  //
  Grid* bias_grid_pntr_;
  Grid* fes_grid_pntr_;
  Grid* ps_grid_pntr_;
  //
 public:
  static void registerKeywords( Keywords& keys );
  // Constructor
  explicit LinearBasisSetExpansion(
    const std::string&,
    Communicator&,
    std::vector<Value*>,
    std::vector<BasisFunctions*>,
    CoeffsVector* bias_coeffs_pntr_in=NULL);
  //
  ~LinearBasisSetExpansion();
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
  CoeffsVector& CoeffDerivsAverTargetDist() const;
  //
  void setSerial() {serial_=true;}
  void setParallel() {serial_=false;}
  //
  void linkVesBias(bias::VesBias*);
  void linkAction(Action*);
  // calculate bias and derivatives
  static double getBiasAndForces(const std::vector<double>&, std::vector<double>&, std::vector<double>&, std::vector<BasisFunctions*>&, CoeffsVector*, Communicator* comm_in=NULL);
  double getBiasAndForces(const std::vector<double>&, std::vector<double>&, std::vector<double>&);
  double getBias(const std::vector<double>&);
  double getFES_WellTempered(const std::vector<double>&);
  //
  static void getBasisSetValues(const std::vector<double>&, std::vector<double>&, std::vector<BasisFunctions*>&, CoeffsVector*, Communicator* comm_in=NULL);
  void getBasisSetValues(const std::vector<double>&, std::vector<double>&);
  // Grid stuff
  void setupBiasGrid(const bool usederiv=false);
  void updateBiasGrid();
  void writeBiasGridToFile(const std::string&, const bool);
  //
  std::vector<unsigned int> getGridBins() const {return grid_bins_;}
  void setGridBins(const std::vector<unsigned int>&);
  void setGridBins(const unsigned int);
  //
  void setupUniformTargetDistribution();
  //
  void setupTargetDistribution(const std::vector<TargetDistribution*>&);
  void setupTargetDistribution(const std::vector<std::string>&);
  void setupTargetDistribution(TargetDistribution*);
  void setupTargetDistribution(const std::string&);
  //
  // Well-Tempered p(s) stuff
  void setupWellTemperedTargetDistribution(const double, const std::vector<unsigned int>&);
  void updateWellTemperedFESCoeffs();
  double getWellTemperedBiasFactor() const {return welltemp_biasf_;}
private:
  //
  Grid* setupGeneralGrid(const std::vector<unsigned int>&, const bool usederiv=false);
  void setupSeperableTargetDistribution(const std::vector<TargetDistribution*>&);
  void setupNonSeperableTargetDistribution(const TargetDistribution*);
  //
  void calculateCoeffDerivsAverFromGrid(const Grid*, const bool normalize_dist=false);
  //
};

inline
std::vector<Value*> LinearBasisSetExpansion::getPntrsToArguments() const {return args_pntrs_;}

inline
std::vector<BasisFunctions*> LinearBasisSetExpansion::getPntrsToBasisFunctions() const {return basisf_pntrs_;}

inline
CoeffsVector* LinearBasisSetExpansion::getPntrToBiasCoeffs() const {return bias_coeffs_pntr_;}

inline
Grid* LinearBasisSetExpansion::getPntrToBiasGrid() const {return bias_grid_pntr_;}

inline
unsigned int LinearBasisSetExpansion::getNumberOfArguments() const {return nargs_;}

inline
std::vector<unsigned int> LinearBasisSetExpansion::getNumberOfBasisFunctions() const {return nbasisf_;}

inline
size_t LinearBasisSetExpansion::getNumberOfCoeffs() const {return ncoeffs_;}

inline
CoeffsVector& LinearBasisSetExpansion::BiasCoeffs() const {return *bias_coeffs_pntr_;}

inline
CoeffsVector& LinearBasisSetExpansion::FesWTCoeffs() const {return *fes_wt_coeffs_pntr_;}

inline
CoeffsVector& LinearBasisSetExpansion::CoeffDerivsAverTargetDist() const {return *coeffderivs_aver_ps_pntr_;}

inline
void LinearBasisSetExpansion::setupTargetDistribution(const std::string& targetdist_keyword) {
  std::vector<std::string> targetdist_keywords(1);
  targetdist_keywords[0] = targetdist_keyword;
  setupTargetDistribution(targetdist_keywords);
}

inline
void LinearBasisSetExpansion::setupTargetDistribution(TargetDistribution* targetdist_pntr) {
  std::vector<TargetDistribution*> targetdist_pntrs(1);
  targetdist_pntrs[0] = targetdist_pntr;
  setupTargetDistribution(targetdist_pntrs);
}


inline
double LinearBasisSetExpansion::getBiasAndForces(const std::vector<double>& args_values, std::vector<double>& forces, std::vector<double>& coeffsderivs_values) {
  return getBiasAndForces(args_values,forces,coeffsderivs_values,basisf_pntrs_, bias_coeffs_pntr_, &mycomm_);
}


inline
double LinearBasisSetExpansion::getBias(const std::vector<double>& args_values) {
  std::vector<double> forces_dummy(nargs_);
  std::vector<double> coeffsderivs_values_dummy(ncoeffs_);
  return getBiasAndForces(args_values,forces_dummy,coeffsderivs_values_dummy,basisf_pntrs_, bias_coeffs_pntr_, &mycomm_);
}


inline
double LinearBasisSetExpansion::getFES_WellTempered(const std::vector<double>& args_values) {
  std::vector<double> forces_dummy(nargs_);
  std::vector<double> coeffsderivs_values_dummy(ncoeffs_);
  return getBiasAndForces(args_values,forces_dummy,coeffsderivs_values_dummy,basisf_pntrs_, fes_wt_coeffs_pntr_, &mycomm_);
}


inline
void LinearBasisSetExpansion::getBasisSetValues(const std::vector<double>& args_values, std::vector<double>& basisset_values) {
  getBasisSetValues(args_values,basisset_values,basisf_pntrs_, bias_coeffs_pntr_, &mycomm_);
}


inline
void LinearBasisSetExpansion::setupBiasGrid(const bool usederiv) {
  bias_grid_pntr_ = setupGeneralGrid(grid_bins_,usederiv);
}


}

}

#endif
