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


#define PLUMED_VESBIAS_INIT(ao) Action(ao),VesBias(ao)

namespace PLMD{

  class CoeffsVector;
  class CoeffsMatrix;
  class BasisFunctions;
  class Value;
  class Optimizer;

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
  CoeffsVector* coeffs_pntr;
  CoeffsVector* coeffderivs_aver_ps_pntr;
  CoeffsVector* gradient_pntr;
  CoeffsMatrix* hessian_pntr;
  std::vector<double> coeffderivs_aver_sampled;
  std::vector<double> coeffderivs_cov_sampled;
  //
  Optimizer* optimizer_pntr;
  bool optimize_coeffs_;
  //
  bool compute_hessian_;
  bool diagonal_hessian_;
  //
  double aver_counter;
  double kbt_;
private:
  void initializeGradientAndHessian();
protected:
  void initializeCoeffs(const std::vector<std::string>&, const std::vector<unsigned int>&);
  void initializeCoeffs(std::vector<Value*>&, std::vector<BasisFunctions*>&);
  void linkCoeffs(CoeffsVector*);
  void linkCoeffs(CoeffsVector&);
  void setCoeffsDerivs(const std::vector<double>&);
  void setCoeffsDerivsOverTargetDist(const std::vector<double>&);
public:
  static void registerKeywords(Keywords&);
  explicit VesBias(const ActionOptions&ao);
  ~VesBias();
  //
  CoeffsVector* getCoeffsPtr() const {return coeffs_pntr;}
  CoeffsVector* getCoeffDerivsAverTargetDistPtr() const {return coeffderivs_aver_ps_pntr;}
  CoeffsVector* getGradientPtr()const {return gradient_pntr;}
  CoeffsMatrix* getHessianPtr() const {return hessian_pntr;}
  //
  size_t numberOfCoeffs() const;
  double getKbT() const;
  double getBeta() const;
  //
  CoeffsVector& Coeffs() const;
  CoeffsVector& CoeffDerivsAverTargetDist() const;
  CoeffsVector& Gradient() const;
  CoeffsMatrix& Hessian() const;
  //
  size_t getCoeffsIndex(const std::vector<unsigned int>& indices) const;
  std::vector<unsigned int> getCoeffsIndices(const size_t index) const;
  size_t getHessianIndex(const size_t index1, const size_t index2) const;
  //
  bool computeHessian() const {return compute_hessian_;}
  bool diagonalHessian() const {return diagonal_hessian_;}
  //
  bool optimizeCoeffs() const {return optimize_coeffs_;}
  //
  void updateGradientAndHessian();
  void clearGradientAndHessian();
  //
  void linkOptimizer(Optimizer*);
  void enableHessian(const bool diagonal_hessian=true);
  void disableHessian();
};

inline
size_t VesBias::numberOfCoeffs() const {return coeffs_pntr->numberOfCoeffs();}

inline
CoeffsVector& VesBias::Coeffs() const {return *coeffs_pntr;}

inline
CoeffsVector& VesBias::CoeffDerivsAverTargetDist() const {return *coeffderivs_aver_ps_pntr;}

inline
CoeffsVector& VesBias::Gradient() const {return *gradient_pntr;}

inline
CoeffsMatrix& VesBias::Hessian() const {return *hessian_pntr;}

inline
double VesBias::getKbT() const {return kbt_;}

inline
double VesBias::getBeta() const {return 1.0/kbt_;}

inline
size_t VesBias::getCoeffsIndex(const std::vector<unsigned int>& indices) const {return coeffs_pntr->getIndex(indices);}

inline
std::vector<unsigned int> VesBias::getCoeffsIndices(const size_t index) const {return coeffs_pntr->getIndices(index);}

inline
size_t VesBias::getHessianIndex(const size_t index1, const size_t index2) const {return hessian_pntr->getMatrixIndex(index1,index2);}


}
}

#endif
