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
  class TargetDistribution;

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
  std::vector<CoeffsVector*> coeffs_pntrs_;
  std::vector<CoeffsVector*> coeffderivs_aver_ps_pntrs_;
  std::vector<CoeffsVector*> gradient_pntrs_;
  std::vector<CoeffsMatrix*> hessian_pntrs_;
  std::vector<std::vector<double> > coeffderivs_aver_sampled;
  std::vector<std::vector<double> >coeffderivs_cov_sampled;
  bool use_multiple_coeffssets_;
  //
  std::vector<std::string> coeffs_fnames;
  //
  size_t ncoeffs_total_;
  //
  Optimizer* optimizer_pntr_;
  bool optimize_coeffs_;
  //
  bool compute_hessian_;
  bool diagonal_hessian_;
  //
  std::vector<std::string> targetdist_keywords_;
  std::vector<TargetDistribution*> targetdist_pntrs_;
  bool dynamic_targetdist_;
  //
  std::string fname_coeffderivs_aver_ps;
  //
  double aver_counter;
  double kbt_;
  //
  bool uniform_targetdist_;
  //
  double welltemp_biasf_;
  bool welltemp_targetdist_;
  //
  std::vector<unsigned int> grid_bins_;
  std::vector<double> grid_min_;
  std::vector<double> grid_max_;
  //
  std::string bias_filename_;
  std::string fes_filename_;
  std::string targetdist_filename_;
  //
private:
  void initializeCoeffs(CoeffsVector*);
protected:
  //
  void checkThatTemperatureIsGiven();
  //
  void addCoeffsSet(const std::vector<std::string>&,const std::vector<unsigned int>&);
  void addCoeffsSet(std::vector<Value*>&,std::vector<BasisFunctions*>&);
  void addCoeffsSet(CoeffsVector*);
  std::string getCoeffsSetLabelString(const std::string&, const unsigned int coeffs_id = 0);
  void clearCoeffsPntrsVector();
  void setCoeffsDerivs(const std::vector<double>&, const unsigned int c_id = 0);
  void setCoeffsDerivsOverTargetDist(const std::vector<double>&, const unsigned int coeffs_id = 0);
  void setCoeffsDerivsOverTargetDist(const CoeffsVector&, const unsigned coeffs_id= 0);
  void setCoeffsDerivsOverTargetDistToZero(const unsigned coeffs_id= 0);

  void readCoeffsFromFiles();
  //
  template<class T>
  bool parseMultipleValues(const std::string&, std::vector<T>&, unsigned int);
  template<class T>
  bool parseMultipleValues(const std::string&, std::vector<T>&, unsigned int, const T&);
public:
  static void registerKeywords(Keywords&);
  explicit VesBias(const ActionOptions&ao);
  ~VesBias();
  //
  void apply();
  //
  static void useInitialCoeffsKeywords(Keywords&);
  static void useTargetDistributionKeywords(Keywords&);
  static void useWellTemperdKeywords(Keywords&);
  static void useGridBinKeywords(Keywords&);
  static void useGridLimitsKeywords(Keywords&);
  //
  std::vector<CoeffsVector*> getCoeffsPntrs() const {return coeffs_pntrs_;}
  std::vector<CoeffsVector*> getCoeffDerivsAverTargetDistPntrs() const {return coeffderivs_aver_ps_pntrs_;}
  std::vector<CoeffsVector*> getGradientPntrs()const {return gradient_pntrs_;}
  std::vector<CoeffsMatrix*> getHessianPntrs() const {return hessian_pntrs_;}
  std::vector<TargetDistribution*> getTargetDistributionPntrs() const {return targetdist_pntrs_;}
  //
  CoeffsVector* getCoeffsPntr(const unsigned int coeffs_id = 0) const {return coeffs_pntrs_[coeffs_id];}
  CoeffsVector* getCoeffDerivsAverTargetDistPntr(const unsigned int coeffs_id = 0) const {return coeffderivs_aver_ps_pntrs_[coeffs_id];}
  CoeffsVector* getGradientPntr(const unsigned int coeffs_id = 0)const {return gradient_pntrs_[coeffs_id];}
  CoeffsMatrix* getHessianPntr(const unsigned int coeffs_id = 0) const {return hessian_pntrs_[coeffs_id];}
  //
  std::vector<std::string> getTargetDistributionKeywords() const {return targetdist_keywords_;}
  unsigned int getNumberOfTargetDistributionKeywords() const {return targetdist_keywords_.size();}
  //
  size_t numberOfCoeffs(const unsigned int coeffs_id = 0) const;
  size_t totalNumberOfCoeffs() const;
  unsigned int numberOfCoeffsSets() const;
  double getKbT() const {return kbt_;}
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
  Optimizer* getOptimizerPntr() const {return optimizer_pntr_;}

  void updateGradientAndHessian();
  void clearGradientAndHessian();
  //
  virtual void updateTargetDistributions();
  void writeCoeffDerivsAverTargetDistToFile(const bool append = false, const unsigned int iteration = 0);
  //
  void linkOptimizer(Optimizer*);
  void enableHessian(const bool diagonal_hessian=true);
  void disableHessian();
  //
  void enableMultipleCoeffsSets() {use_multiple_coeffssets_=true;}
  //
  void enableDynamicTargetDistribution() {dynamic_targetdist_=true;}
  void disableDynamicTargetDistribution() {dynamic_targetdist_=false;}
  bool dynamicTargetDistribution() const {return dynamic_targetdist_;}
  //
  void enableUniformTargetDistribution() {uniform_targetdist_=true;}
  void disableUniformTargetDistribution() {uniform_targetdist_=false;}
  bool uniformTargetDistribution() const {return uniform_targetdist_;}
    //
  double getWellTemperedBiasFactor() const;
  bool wellTemperdTargetDistribution() const {return welltemp_targetdist_;}
  //
  std::vector<unsigned int> getGridBins() const {return grid_bins_;}
  void setGridBins(const std::vector<unsigned int>&);
  void setGridBins(const unsigned int);
  std::vector<double> getGridMax() const {return grid_max_;}
  void setGridMax(const std::vector<double>&);
  std::vector<double> getGridMin() const {return grid_min_;}
  void setGridMin(const std::vector<double>&);
  //
  std::string getBiasOutputFilename() const {return bias_filename_;}
  std::string getCurrentBiasOutputFilename() const;
  std::string getFesOutputFilename() const {return fes_filename_;}
  std::string getCurrentFesOutputFilename() const;
  std::string getTargetDistOutputFilename() const {return targetdist_filename_;}
  std::string getCurrentTargetDistOutputFilename(const std::string suffix="") const;
  //
  virtual void setupBiasFileOutput() {};
  virtual void writeBiasToFile() {};
  virtual void setupFesFileOutput() {};
  virtual void writeFesToFile() {};
};

inline
size_t VesBias::numberOfCoeffs(const unsigned int coeffs_id) const {return coeffs_pntrs_[coeffs_id]->numberOfCoeffs();}

inline
unsigned int VesBias::numberOfCoeffsSets() const {return ncoeffssets_;}

inline
size_t VesBias::totalNumberOfCoeffs() const {return ncoeffs_total_;}

inline
CoeffsVector& VesBias::Coeffs(const unsigned int coeffs_id) const {return *coeffs_pntrs_[coeffs_id];}

inline
CoeffsVector& VesBias::CoeffDerivsAverTargetDist(const unsigned int coeffs_id) const {return *coeffderivs_aver_ps_pntrs_[coeffs_id];}

inline
CoeffsVector& VesBias::Gradient(const unsigned int coeffs_id) const {return *gradient_pntrs_[coeffs_id];}

inline
CoeffsMatrix& VesBias::Hessian(const unsigned int coeffs_id) const {return *hessian_pntrs_[coeffs_id];}

inline
size_t VesBias::getCoeffsIndex(const std::vector<unsigned int>& indices, const unsigned int coeffs_id) const {return coeffs_pntrs_[coeffs_id]->getIndex(indices);}

inline
std::vector<unsigned int> VesBias::getCoeffsIndices(const size_t index, const unsigned int coeffs_id) const {return coeffs_pntrs_[coeffs_id]->getIndices(index);}

inline
size_t VesBias::getHessianIndex(const size_t index1, const size_t index2, const unsigned int coeffs_id) const {return hessian_pntrs_[coeffs_id]->getMatrixIndex(index1,index2);}

inline
double VesBias::getWellTemperedBiasFactor() const {
  plumed_massert(welltemp_targetdist_,"the well-tempered target distribution is not active so it doesn't make sense to get the value of the bias factor");
  return welltemp_biasf_;
}

inline
double VesBias::getBeta() const {
  plumed_massert(kbt_!=0.0,"you are requesting beta=1/(kB*T) when kB*T has not been defined. You need to give the temperature using the TEMP keyword as the MD engine does not pass it to PLUMED.");
  return 1.0/kbt_;
}


template<class T>
bool VesBias::parseMultipleValues(const std::string& keyword, std::vector<T>& values, unsigned int nvalues) {
  plumed_assert(nvalues>0);
  plumed_assert(values.size()==0);
  bool identical_values=false;
  //
  parseVector(keyword,values);
  if(values.size()==1 && nvalues>1){
    values.resize(nvalues,values[0]);
    identical_values=true;
  }
  if(values.size()>0 && values.size()!=nvalues){
    std::string s1; Tools::convert(nvalues,s1);
    plumed_merror("Error in " + keyword + " keyword: either give 1 common parameter value or " + s1 + " seperate parameter values");
  }
  return identical_values;
}

template<class T>
bool VesBias::parseMultipleValues(const std::string& keyword, std::vector<T>& values, unsigned int nvalues, const T& default_value) {
  bool identical_values = parseMultipleValues(keyword,values,nvalues);
  if(values.size()==0){
    values.resize(nvalues,default_value);
    identical_values=true;
  }
  return identical_values;
}



}
}

#endif
