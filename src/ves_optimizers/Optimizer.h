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

class Optimizer :
 public ActionPilot,
 public ActionWithValue
{
private:
  bool usehessian_;
  std::string description_;
  std::string type_;
  double step_size_;
  double current_step_size_;
  //
  CoeffsVector* coeffs_ptr;
  CoeffsVector* aux_coeffs_ptr;
  CoeffsVector* gradient_ptr;
  CoeffsMatrix* hessian_ptr;
  //
  bias::VesBias* bias_ptr;
  //
private:
  void updateOutputComponents();
protected:
  void turnOnHessian();
  void turnOffHessian();
  CoeffsVector& Coeffs() const;
  CoeffsVector& AuxCoeffs() const;
  CoeffsVector& Gradient() const;
  CoeffsMatrix& Hessian() const;
  double StepSize() const;
  virtual void coeffsUpdate()=0;
  void setCurrentStepSize(const double);
public:
  static void registerKeywords(Keywords&);
  Optimizer(const ActionOptions&ao);
  ~Optimizer();
  std::string getType() const;
  std::string getDescription() const;
  //
  double getStepSize() const;
  double getCurrentStepSize() const;
  void setStepSize(const double step_size_in){step_size_ = step_size_in;}
  //
  void apply(){};
  void calculate(){};
  void update();
  unsigned int getNumberOfDerivatives(){return 0;}
  //
  bool useHessian() const {return usehessian_;}  ;
  };

inline
std::string Optimizer::getType() const {return type_;}

inline
std::string Optimizer::getDescription() const {return description_;}

inline
double Optimizer::StepSize() const {return step_size_;}

inline
CoeffsVector& Optimizer::Coeffs() const {return *coeffs_ptr;}

inline
CoeffsVector& Optimizer::AuxCoeffs() const {return *aux_coeffs_ptr;}

inline
CoeffsVector& Optimizer::Gradient() const {return *gradient_ptr;}

inline
CoeffsMatrix& Optimizer::Hessian() const {return *hessian_ptr;}

inline
double Optimizer::getStepSize() const {return step_size_;}

inline
double Optimizer::getCurrentStepSize() const {return current_step_size_;}


}

#endif