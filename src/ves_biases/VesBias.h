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
  unsigned int ncoeffssets_;
  std::vector<CoeffsVector*> coeffs_pntrs;
  std::vector<CoeffsVector*> coeffderivs_aver_ps_pntrs;
  std::vector<CoeffsVector*> gradient_pntrs;
  std::vector<CoeffsMatrix*> hessian_pntrs;
  std::vector<std::vector<double> > coeffderivs_aver_sampled;
  std::vector<std::vector<double> >coeffderivs_cov_sampled;
  bool use_multiple_coeffssets_;
  //
  std::vector<std::string> coeffs_fnames;
  //
  size_t ncoeffs_total_;
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
  void initializeCoeffs(CoeffsVector*);
  std::string labelString(const std::string&, const unsigned int coeffs_id = 0);
protected:
  void addCoeffsSet(const std::vector<std::string>&,const std::vector<unsigned int>&);
  void addCoeffsSet(std::vector<Value*>&,std::vector<BasisFunctions*>&);
  void setCoeffsDerivs(const std::vector<double>&, const unsigned int c_id = 0);
  void setCoeffsDerivsOverTargetDist(const std::vector<double>&, const unsigned int coeffs_id = 0);
  void readCoeffsFromFiles();
public:
  static void registerKeywords(Keywords&);
  explicit VesBias(const ActionOptions&ao);
  ~VesBias();
  //
  void apply();
  //
  std::vector<CoeffsVector*> getCoeffsPntrs() const {return coeffs_pntrs;}
  std::vector<CoeffsVector*> getCoeffDerivsAverTargetDistPntrs() const {return coeffderivs_aver_ps_pntrs;}
  std::vector<CoeffsVector*> getGradientPntrs()const {return gradient_pntrs;}
  std::vector<CoeffsMatrix*> getHessianPntrs() const {return hessian_pntrs;}
  //
  CoeffsVector* getCoeffsPntr(const unsigned int coeffs_id = 0) const {return coeffs_pntrs[coeffs_id];}
  CoeffsVector* getCoeffDerivsAverTargetDistPntr(const unsigned int coeffs_id = 0) const {return coeffderivs_aver_ps_pntrs[coeffs_id];}
  CoeffsVector* getGradientPntr(const unsigned int coeffs_id = 0)const {return gradient_pntrs[coeffs_id];}
  CoeffsMatrix* getHessianPntr(const unsigned int coeffs_id = 0) const {return hessian_pntrs[coeffs_id];}

  //
  size_t numberOfCoeffs(const unsigned int coeffs_id = 0) const;
  size_t totalNumberOfCoeffs() const;
  unsigned int numberOfCoeffsSets() const;
  double getKbT() const;
  double getBeta() const;
  //
  CoeffsVector& Coeffs(const unsigned int coeffs_id = 0) const;
  CoeffsVector& CoeffDerivsAverTargetDist(const unsigned int coeffs_id = 0) const;
  CoeffsVector& Gradient(const unsigned int coeffs_id = 0) const;
  CoeffsMatrix& Hessian(const unsigned int coeffs_id = 0) const;
  //
  size_t getCoeffsIndex(const std::vector<unsigned int>& indices, const unsigned int coeffs_id = 0) const;
  std::vector<unsigned int> getCoeffsIndices(const size_t index, const unsigned int coeffs_id = 0) const;
  size_t getHessianIndex(const size_t index1, const size_t index2, const unsigned int coeffs_id = 0) const;
  //
  bool computeHessian() const {return compute_hessian_;}
  bool diagonalHessian() const {return diagonal_hessian_;}
  //
  bool optimizeCoeffs() const {return optimize_coeffs_;}

  void updateGradientAndHessian();
  void clearGradientAndHessian();
  //
  void linkOptimizer(Optimizer*);
  void enableHessian(const bool diagonal_hessian=true);
  void disableHessian();
  //
  void enableMultipleCoeffsSets() {use_multiple_coeffssets_=true;}
};

inline
size_t VesBias::numberOfCoeffs(const unsigned int coeffs_id) const {return coeffs_pntrs[coeffs_id]->numberOfCoeffs();}

inline
unsigned int VesBias::numberOfCoeffsSets() const {return ncoeffssets_;}

inline
size_t VesBias::totalNumberOfCoeffs() const {return ncoeffs_total_;}

inline
CoeffsVector& VesBias::Coeffs(const unsigned int coeffs_id) const {return *coeffs_pntrs[coeffs_id];}

inline
CoeffsVector& VesBias::CoeffDerivsAverTargetDist(const unsigned int coeffs_id) const {return *coeffderivs_aver_ps_pntrs[coeffs_id];}

inline
CoeffsVector& VesBias::Gradient(const unsigned int coeffs_id) const {return *gradient_pntrs[coeffs_id];}

inline
CoeffsMatrix& VesBias::Hessian(const unsigned int coeffs_id) const {return *hessian_pntrs[coeffs_id];}

inline
double VesBias::getKbT() const {return kbt_;}

inline
double VesBias::getBeta() const {return 1.0/kbt_;}

inline
size_t VesBias::getCoeffsIndex(const std::vector<unsigned int>& indices, const unsigned int coeffs_id) const {return coeffs_pntrs[coeffs_id]->getIndex(indices);}

inline
std::vector<unsigned int> VesBias::getCoeffsIndices(const size_t index, const unsigned int coeffs_id) const {return coeffs_pntrs[coeffs_id]->getIndices(index);}

inline
size_t VesBias::getHessianIndex(const size_t index1, const size_t index2, const unsigned int coeffs_id) const {return hessian_pntrs[coeffs_id]->getMatrixIndex(index1,index2);}


}
}

#endif
