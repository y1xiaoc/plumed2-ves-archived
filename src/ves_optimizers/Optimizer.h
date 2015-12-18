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
  std::string description_;
  std::string type_;
  //
  double stepsize_;
  double current_stepsize_;
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
  unsigned int coeffs_wstride_;
  std::string coeffs_fname_;
  OFile coeffsOfile_;
  //
  unsigned int gradient_wstride_;
  std::string gradient_fname_;
  OFile gradientOfile_;
  //
  unsigned int hessian_wstride_;
  std::string hessian_fname_;
  OFile hessianOfile_;
  //
  CoeffsVector* coeffs_ptr;
  CoeffsVector* aux_coeffs_ptr;
  CoeffsVector* gradient_ptr;
  CoeffsMatrix* hessian_ptr;
  CoeffsVector* coeffs_mask_ptr;
  //
  bias::VesBias* bias_ptr;
  //
private:
  void updateOutputComponents();
  void writeOutputFiles();
protected:
  void turnOnHessian();
  void turnOffHessian();
  CoeffsMatrix* switchToDiagonalHessian(bias::VesBias*);
  CoeffsMatrix* switchToFullHessian(bias::VesBias*);
  //
  CoeffsVector& Coeffs() const;
  CoeffsVector& AuxCoeffs() const;
  CoeffsVector& Gradient() const;
  CoeffsMatrix& Hessian() const;
  CoeffsVector& CoeffsMask() const;
  double StepSize() const;
  virtual void coeffsUpdate()=0;
  void setCurrentStepSize(const double);

public:
  static void registerKeywords(Keywords&);
  static void useMultipleWalkersKeywords(Keywords&);
  static void useHessianKeywords(Keywords&);
  static void useFixedStepSizeKeywords(Keywords&);
  static void useChangingStepSizeKeywords(Keywords&);
  static void useMaskKeywords(Keywords&);
  //
  explicit Optimizer(const ActionOptions&ao);
  ~Optimizer();
  std::string getType() const;
  std::string getDescription() const;
  //
  double getStepSize() const;
  double getCurrentStepSize() const;
  void setStepSize(const double stepsize_in){stepsize_ = stepsize_in;}
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
  //
  bool useHessian() const {return use_hessian_;}
  bool diagonalHessian() const {return diagonal_hessian_;}
  //
  bool useMultipleWalkers() const {return use_mwalkers_mpi_;}
  //
  CoeffsVector* getCoeffsPtr() const {return coeffs_ptr;}
  CoeffsVector* getAuxCoeffsPtr() const {return aux_coeffs_ptr;}
  CoeffsVector* getGradientPtr()const {return gradient_ptr;}
  CoeffsMatrix* getHessianPtr() const {return hessian_ptr;}
  CoeffsVector* getCoeffsMaskPtr() const {return coeffs_mask_ptr;}
  };

inline
std::string Optimizer::getType() const {return type_;}

inline
std::string Optimizer::getDescription() const {return description_;}

inline
double Optimizer::StepSize() const {return stepsize_;}

inline
CoeffsVector& Optimizer::Coeffs() const {return *coeffs_ptr;}

inline
CoeffsVector& Optimizer::AuxCoeffs() const {return *aux_coeffs_ptr;}

inline
CoeffsVector& Optimizer::Gradient() const {return *gradient_ptr;}

inline
CoeffsMatrix& Optimizer::Hessian() const {
  plumed_massert(use_hessian_,"You cannot use the Hessian without asking for before");
  return *hessian_ptr;
}

inline
CoeffsVector& Optimizer::CoeffsMask() const {return *coeffs_mask_ptr;}

inline
double Optimizer::getStepSize() const {return stepsize_;}

inline
double Optimizer::getCurrentStepSize() const {return current_stepsize_;}

inline
unsigned int Optimizer::getIterationCounter() const {return iter_counter;}

inline
double Optimizer::getIterationCounterDbl() const {return (double) iter_counter;}

inline
void Optimizer::increaseIterationCounter() {iter_counter++;}


}

#endif
