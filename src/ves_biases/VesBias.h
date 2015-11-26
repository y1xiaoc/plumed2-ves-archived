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
#ifndef __PLUMED_ves_biases_VesBias_h
#define __PLUMED_ves_biases_VesBias_h

#include "ves_tools/CoeffsVector.h"
#include "ves_tools/CoeffsMatrix.h"

#include "core/ActionPilot.h"
#include "core/ActionWithValue.h"
#include "core/ActionWithArguments.h"
#include "bias/Bias.h"

#include <vector>
#include <string>
#include <cmath>


#define PLUMED_VARIATIONALBIAS_INIT(ao) Action(ao),VesBias(ao)

namespace PLMD{

  class CoeffsVector;
  class CoeffsMatrix;
  class BasisFunctions;
  class Value;

namespace bias{

/**
\ingroup INHERIT
Abstract base class for implementing biases the extents the normal Bias.h class
to include functions related to the variational approach.
*/

class VesBias:
public Bias
{
private:
  CoeffsVector* coeffs_ptr;
  CoeffsVector* coeffderivs_aver_ps_ptr;
  CoeffsVector* gradient_ptr;
  CoeffsMatrix* hessian_ptr;
  std::vector<double> coeffderivs_aver_sampled;
  std::vector<double> coeffderivs_cov_sampled;
  //
  bool hessian_diagonal_;
  bool use_mwalkers_;
  //
  double aver_counter;
private:
  void initializeGradientAndHessian();
protected:
  void initializeCoeffs(const std::vector<std::string>&, const std::vector<unsigned int>&);
  void initializeCoeffs(std::vector<Value*>&, std::vector<BasisFunctions*>&);
  void linkCoeffs(CoeffsVector*);
  void linkCoeffs(CoeffsVector&);
  void setCoeffsDerivs(const std::vector<double>&);
public:
  static void registerKeywords(Keywords&);
  VesBias(const ActionOptions&ao);
  ~VesBias();
  //
  CoeffsVector* getCoeffsPtr() const {return coeffs_ptr;}
  CoeffsVector* getCoeffDerivsAverTargetDistPtr() const {return coeffderivs_aver_ps_ptr;}
  CoeffsVector* getGradientPtr()const {return gradient_ptr;}
  CoeffsMatrix* getHessianPtr() const {return hessian_ptr;}
  //
  size_t numberOfCoeffs() const;
  //
  CoeffsVector& Coeffs() const;
  CoeffsVector& CoeffDerivsAverTargetDist() const;
  CoeffsVector& Gradient() const;
  CoeffsMatrix& Hessian() const;
  //
  bool diagonalHessian() const {return hessian_diagonal_;}
  bool useMultipleWalkers() const {return use_mwalkers_;}
  //
  void updateGradientAndHessian();
  void clearGradientAndHessian();
};

inline
size_t VesBias::numberOfCoeffs() const {return coeffs_ptr->numberOfCoeffs();}

inline
CoeffsVector& VesBias::Coeffs() const {return *coeffs_ptr;}

inline
CoeffsVector& VesBias::CoeffDerivsAverTargetDist() const {return *coeffderivs_aver_ps_ptr;}

inline
CoeffsVector& VesBias::Gradient() const {return *gradient_ptr;}

inline
CoeffsMatrix& VesBias::Hessian() const {return *hessian_ptr;}

}
}

#endif
