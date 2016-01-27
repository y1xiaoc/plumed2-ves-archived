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
#ifndef __PLUMED_ves_basisfunctions_BasisFunctions_h
#define __PLUMED_ves_basisfunctions_BasisFunctions_h

#include <vector>
#include <string>
#include <cmath>
#include "core/ActionWithValue.h"

#define PLUMED_BASISFUNCTIONS_INIT(ao) Action(ao),BasisFunctions(ao)

namespace PLMD{

/**
\ingroup INHERIT
Abstract base class for implenting new 1D basis sets.
*/

class TargetDistribution;

class BasisFunctions :
  public ActionWithValue
{
private:
  // print extra info about the basis set
  bool print_debug_info_;
  // to check if the basis set has been defined
  bool has_been_set;
  // description of the basis set
  std::string description_;
  // the type of the basis set
  std::string type_;
  // the maximum order of the basis functions
  unsigned int norder_;
  // the total number of basis functions
  unsigned int nbasis_;
  // the keywords used to invoke the basis set
  std::vector<std::string> bf_keywords_;
  // prefix for the basis function labels
  std::string bf_label_prefix_;
  // label of each basis function
  std::vector<std::string> bf_labels_;
  // if the basis functions are periodic or not
  bool periodic_;
  // if the basis functions are defined on a bounded interval or not
  bool interval_bounded_;
  // the intrinsic interval of the basis functions
  double interval_intrinsic_min_;
  double interval_intrinsic_max_;
  double interval_intrinsic_range_;
  double interval_intrinsic_mean_;
  // the defined (translated) interval of the basis functions
  double interval_min_;
  double interval_max_;
  double interval_range_;
  double interval_mean_;
  // the derivative term in the chain rule coming from the translation of the interval
  double argT_derivf_;
  // calculate numerically the integrals of the basis functions over the intervals
  bool numerical_uniform_integrals_;
  unsigned int nbins_;
  // the integrals of the basis functions over the interval on which they are defined
  std::vector <double> uniform_integrals_;
protected:
  // setup various stuff
  void setupBF();
  void setupInterval();
  void setNumericalIntegrationBins(const unsigned int nbins) {nbins_=nbins;}
  void numericalUniformIntegrals();
  std::vector<double> numericalTargetDistributionIntegrals(const TargetDistribution*) const ;
  virtual void setupLabels();
  virtual void setupUniformIntegrals();
  template<typename T>
  void addKeywordToList(const std::string&, const T);
  void addKeywordToList(const std::string&, const bool);
  //
  void setPeriodic() {periodic_=true;}
  void setNonPeriodic() {periodic_=false;}
  void setIntervalBounded() {interval_bounded_=true;}
  void setIntervalNonBounded() {interval_bounded_=false;}
  void setType(const std::string& type_in) {type_=type_in;}
  void setDescription(const std::string& description_in) {description_=description_in;}
  //
  void setNumberOfBasisFunctions(const unsigned int);
  void setIntrinsicInterval(const double, const double);
  //
  void setUniformIntegral(const unsigned int, const double);
  void setUniformIntegrals(const std::vector<double>&);
  void setAllUniformIntegralsToZero();
  //
  void setLabelPrefix(const std::string&);
  void setLabel(const unsigned int, const std::string&);
  void setLabels(const std::vector<std::string>&);

public:
  static void registerKeywords(Keywords&);
  explicit BasisFunctions(const ActionOptions&ao);
  bool hasBeenSet() const;
  std::string getType() const;
  std::string getDescription() const;
  unsigned int getOrder() const;
  unsigned int getNumberOfBasisFunctions() const;
  unsigned int numberOfBasisFunctions() const;
  unsigned int getSize() const;
  bool arePeriodic() const;
  bool intervalBounded() const;
  double intervalMin() const;
  double intervalMax() const;
  double intervalRange() const;
  double intervalMean() const;
  double intervalDerivf() const;
  std::vector<double> getUniformIntegrals() const;
  std::vector<double> getTargetDistributionIntegrals(TargetDistribution*) const;
  unsigned getNumberOfDerivatives(){return 0;}
  std::vector<std::string> getKeywordList() const;
  std::string getBasisFunctionLabel(const unsigned int) const;
  std::vector<std::string> getBasisFunctionLabels() const;
  //
  double translateArgument(const double, bool&) const;
  void apply(){};
  void calculate(){};
  // calculate the value for the n-th basis function
  virtual double getValue(const double, const unsigned int, double&, bool&) const = 0;
  // calcuate the values for all basis functions
  virtual void getAllValues(const double, double&, bool&, std::vector<double>&, std::vector<double>&) const = 0;
  //virtual void get2ndDerivaties(const double, std::vector<double>&)=0;
  void printInfo() const;
  std::string getKeywordString() const;
};


inline
bool BasisFunctions::hasBeenSet() const {return has_been_set;}


inline
std::string BasisFunctions::getType() const {return type_;}


inline
std::string BasisFunctions::getDescription() const {return description_;}


inline
unsigned int BasisFunctions::getOrder() const {return norder_;}


inline
unsigned int BasisFunctions::getNumberOfBasisFunctions() const  {return nbasis_;}


inline
unsigned int BasisFunctions::numberOfBasisFunctions() const  {return nbasis_;}


inline
unsigned int BasisFunctions::getSize() const {return getNumberOfBasisFunctions();}


inline
bool BasisFunctions::arePeriodic() const {return periodic_;}


inline
bool BasisFunctions::intervalBounded() const {return interval_bounded_;}


inline
double BasisFunctions::intervalMin() const {return interval_min_;}


inline
double BasisFunctions::intervalMax() const {return interval_max_;}


inline
double BasisFunctions::intervalRange() const {return interval_range_;}


inline
double BasisFunctions::intervalMean() const {return interval_mean_;}


inline
double BasisFunctions::intervalDerivf() const {return argT_derivf_;}


inline
std::vector<double> BasisFunctions::getUniformIntegrals() const {return uniform_integrals_;}


inline
std::vector<double> BasisFunctions::getTargetDistributionIntegrals(TargetDistribution* targetdist_in) const {
  return numericalTargetDistributionIntegrals(targetdist_in);
}


inline
std::vector<std::string> BasisFunctions::getKeywordList() const {return bf_keywords_;}


inline
std::string BasisFunctions::getBasisFunctionLabel(const unsigned int index) const {return bf_labels_[index];}


inline
std::vector<std::string> BasisFunctions::getBasisFunctionLabels() const {return bf_labels_;}


inline
void BasisFunctions::setUniformIntegral(const unsigned index, const double value) {
  uniform_integrals_[index] = value;
}


inline
void BasisFunctions::setUniformIntegrals(const std::vector<double>& uniform_integrals_in) {
  plumed_assert(uniform_integrals_in.size()==nbasis_);
  uniform_integrals_ = uniform_integrals_in;
}


inline
void BasisFunctions::setAllUniformIntegralsToZero() {
  uniform_integrals_.assign(nbasis_,0.0);
}

inline
void BasisFunctions::setLabelPrefix(const std::string& bf_label_prefix_in) {
  bf_label_prefix_ = bf_label_prefix_in;
}


inline
void BasisFunctions::setLabel(const unsigned int index, const std::string& label) {
  bf_labels_[index] = label;
}


inline
void BasisFunctions::setLabels(const std::vector<std::string>& bf_labels_in) {
  bf_labels_ = bf_labels_in;
}


inline
void BasisFunctions::setIntrinsicInterval(const double interval_intrinsic_min_in, const double interval_intrinsic_max_in) {
  interval_intrinsic_min_ = interval_intrinsic_min_in;
  interval_intrinsic_max_ = interval_intrinsic_max_in;
}


inline
double BasisFunctions::translateArgument(const double arg, bool& inside_interval) const {
  inside_interval=true;
  double argT = (arg-interval_mean_)*argT_derivf_;
  if(argT < interval_intrinsic_min_){
    inside_interval=false;
    argT=interval_intrinsic_min_;
  }
  else if(argT > interval_intrinsic_max_){
    inside_interval=false;
    argT=interval_intrinsic_max_;
  }
  return argT;
}


template<typename T>
void BasisFunctions::addKeywordToList(const std::string& keyword, const T value){
  std::string str_value;
  Tools::convert(value,str_value);
  bf_keywords_.push_back(keyword+"="+str_value);
}

inline
void BasisFunctions::addKeywordToList(const std::string& keyword, const bool value){
  if(value){bf_keywords_.push_back(keyword);}
}


}

#endif
