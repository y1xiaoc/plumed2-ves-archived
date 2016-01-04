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
#ifndef __PLUMED_ves_optimizers_Optimizer_h
#define __PLUMED_ves_optimizers_Optimizer_h

#include <vector>
#include <string>
#include <cmath>
#include "core/ActionPilot.h"
#include "core/ActionWithValue.h"
#include "ves_biases/VesBias.h"

#define PLUMED_OPTIMIZER_INIT(ao) Action(ao),Optimizer(ao)

namespace PLMD{

/**
\ingroup INHERIT
Abstract base class for implenting new optimization methods
*/

class CoeffsVector;
class VesBias;
class OFile;

class Optimizer :
 public ActionPilot,
 public ActionWithValue
{
private:
  const std::string description_;
  const std::string type_;
  //
  std::vector<double> stepsizes_;
  std::vector<double> current_stepsizes;
  bool fixed_stepsize_;
  //
  unsigned int iter_counter;
  //
  bool use_hessian_;
  bool diagonal_hessian_;
  //
  bool use_mwalkers_mpi_;
  bool mwalkers_mpi_single_files_;
  //
  std::string fname_prefix_;
  //
  unsigned int coeffs_wstride_;
  std::vector<OFile*> coeffsOFiles_;
  //
  unsigned int gradient_wstride_;
  std::vector<OFile*> gradientOFiles_;
  //
  unsigned int hessian_wstride_;
  std::vector<OFile*> hessianOFiles_;
  //
  unsigned int nbiases_;
  std::vector<bias::VesBias*> bias_pntrs;
  //
  unsigned int ncoeffssets_;
  std::vector<CoeffsVector*> coeffs_pntrs;
  std::vector<CoeffsVector*> aux_coeffs_pntrs;
  std::vector<CoeffsVector*> gradient_pntrs;
  std::vector<CoeffsMatrix*> hessian_pntrs;
  std::vector<CoeffsVector*> coeffs_mask_pntrs;
  //
  bool identical_coeffs_shape_;
  //
private:
  void updateOutputComponents();
  void writeOutputFiles();
  void writeOutputFiles(const unsigned int);
  void setupOFiles(std::vector<std::string>&, std::vector<OFile*>&);
  void readCoeffsFromFiles(const std::vector<std::string>&);
protected:
  void turnOnHessian();
  void turnOffHessian();
  std::vector<CoeffsMatrix*> enableHessian(bias::VesBias*, const bool diagonal_hessian=false);
  // CoeffsMatrix* switchToDiagonalHessian(bias::VesBias*);
  // CoeffsMatrix* switchToFullHessian(bias::VesBias*);
  //
  CoeffsVector& Coeffs(const unsigned int coeffs_id = 0) const;
  CoeffsVector& AuxCoeffs(const unsigned int coeffs_id = 0) const;
  CoeffsVector& Gradient(const unsigned int coeffs_id = 0) const;
  CoeffsMatrix& Hessian(const unsigned int coeffs_id = 0) const;
  CoeffsVector& CoeffsMask(const unsigned int coeffs_id = 0) const;
  double StepSize(const unsigned int coeffs_id = 0) const;
  virtual void coeffsUpdate(const unsigned int coeffs_id = 0) = 0;
  void setCurrentStepSize(const double,const unsigned int i = 0);
  void setCurrentStepSizes(const std::vector<double>);
  //
  void turnOffCoeffsOutputFiles();
  //
  void parseFilenames(const std::string&, std::vector<std::string>&, const std::string& default_fname="");
  template<class T>
  void parseValues(const std::string& keyword, std::vector<T>& values);
public:
  static void registerKeywords(Keywords&);
  static void useMultipleWalkersKeywords(Keywords&);
  static void useHessianKeywords(Keywords&);
  static void useFixedStepSizeKeywords(Keywords&);
  static void useDynamicStepSizeKeywords(Keywords&);
  static void useMaskKeywords(Keywords&);
  //
  explicit Optimizer(const ActionOptions&ao);
  ~Optimizer();
  std::string getType() const {return type_;}
  std::string getDescription() const {return description_;}
  //
  unsigned int numberOfBiases() const {return nbiases_;}
  unsigned int numberOfCoeffsSets() const {return ncoeffssets_;}
  //
  std::vector<double> getStepSizes() const;
  std::vector<double> getCurrentStepSizes() const;
  double getStepSize(const unsigned int coeffs_id = 0) const;
  double getCurrentStepSize(const unsigned int coeffs_id = 0) const;
  void setStepSizes(const std::vector<double>);
  void setStepSize(const double, const unsigned int coeffs_id = 0);
  //
  unsigned int getIterationCounter() const;
  double getIterationCounterDbl() const;
  void setIterationCounter(const unsigned int);
  void increaseIterationCounter();
  //
  void apply(){};
  void calculate(){};
  void update();
  unsigned int getNumberOfDerivatives(){return 0;}
  //
  bool fixedStepSize() const {return fixed_stepsize_;}
  bool dynamicStepSize() const {return !fixed_stepsize_;}
  //
  bool useHessian() const {return use_hessian_;}
  bool diagonalHessian() const {return diagonal_hessian_;}
  //
  bool useMultipleWalkers() const {return use_mwalkers_mpi_;}
  //
  std::vector<CoeffsVector*> getCoeffsPntrs() const {return coeffs_pntrs;}
  std::vector<CoeffsVector*> getAuxCoeffsPntrs() const {return aux_coeffs_pntrs;}
  std::vector<CoeffsVector*> getGradientPntrs()const {return gradient_pntrs;}
  std::vector<CoeffsMatrix*> getHessianPntrs() const {return hessian_pntrs;}
  std::vector<CoeffsVector*> getCoeffsMaskPntrs() const {return coeffs_mask_pntrs;}
  };

inline
double Optimizer::StepSize(const unsigned int coeffs_id) const {return stepsizes_[coeffs_id];}

inline
CoeffsVector& Optimizer::Coeffs(const unsigned int coeffs_id) const {return *coeffs_pntrs[coeffs_id];}

inline
CoeffsVector& Optimizer::AuxCoeffs(const unsigned int coeffs_id) const {return *aux_coeffs_pntrs[coeffs_id];}

inline
CoeffsVector& Optimizer::Gradient(const unsigned int coeffs_id) const {return *gradient_pntrs[coeffs_id];}

inline
CoeffsMatrix& Optimizer::Hessian(const unsigned int coeffs_id) const {
  plumed_massert(use_hessian_,"You cannot use the Hessian without asking for before");
  return *hessian_pntrs[coeffs_id];
}

inline
CoeffsVector& Optimizer::CoeffsMask(const unsigned int coeffs_id) const {return *coeffs_mask_pntrs[coeffs_id];}

inline
std::vector<double> Optimizer::getStepSizes() const {return stepsizes_;}

inline
std::vector<double> Optimizer::getCurrentStepSizes() const {return current_stepsizes;}

inline
double Optimizer::getStepSize(const unsigned int coeffs_id) const {return stepsizes_[coeffs_id];}

inline
double Optimizer::getCurrentStepSize(const unsigned int coeffs_id) const {return current_stepsizes[coeffs_id];}

inline
void Optimizer::setStepSizes(const std::vector<double> stepsizes_in) {
  plumed_assert(stepsizes_in.size()==ncoeffssets_);
  stepsizes_ = stepsizes_in;
}

inline
void Optimizer::setStepSize(const double stepsize_in, const unsigned int coeffs_id) {
  stepsizes_[coeffs_id] = stepsize_in;
}

inline
void Optimizer::setCurrentStepSize(const double current_stepsize_in, const unsigned int coeffs_id) {
  current_stepsizes[coeffs_id] = current_stepsize_in;
}

inline
void Optimizer::setCurrentStepSizes(const std::vector<double> current_stepsizes_in) {
  plumed_assert(current_stepsizes_in.size()==ncoeffssets_);
  current_stepsizes = current_stepsizes_in;
}

inline
unsigned int Optimizer::getIterationCounter() const {return iter_counter;}

inline
double Optimizer::getIterationCounterDbl() const {return static_cast<double>(iter_counter);}

inline
void Optimizer::increaseIterationCounter() {iter_counter++;}

inline
void Optimizer::setIterationCounter(const unsigned int iter_counter_in) {iter_counter = iter_counter_in;}


template<class T>
void Optimizer::parseValues(const std::string& keyword, std::vector<T>& values) {
  plumed_assert(ncoeffssets_>0);
  parseVector(keyword,values);
  if(values.size()==1){
    values.resize(ncoeffssets_,values[0]);
  }
  plumed_massert(values.size()==ncoeffssets_,"Error in " + keyword + " keyword: either give one common value for all coefficient sets or a seperate value for each set");
}

}

#endif
