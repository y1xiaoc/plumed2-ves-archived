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
#ifndef __PLUMED_variational_BasisFunctions_h
#define __PLUMED_variational_BasisFunctions_h

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

class BasisFunctions :
 public ActionWithValue
{
private:

protected:
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
  // prefix for the basis function descriptions 
  std::string bf_description_prefix_; 
  // description of each basis function
  std::vector<std::string> bf_description_; 
  // if the basis functions are periodic or not
  bool periodic_;
  // if the basis functions are defined on a bounded interval or not
  bool interval_bounded_;
  // the default interval of the basis functions 
  double interval_default_min_;
  double interval_default_max_;
  double interval_default_range_;
  double interval_default_mean_;
  // the defined (translated) interval of the basis functions
  double interval_min_;
  double interval_max_;
  double interval_range_;
  double interval_mean_;
  // the derivative term in the chain rule coming from the translation of the interval
  double argT_derivf_;
  // calculate numerically the integrals of the basis functions over the intervals
  bool numerical_bf_integrals_;
  // the integrals of the basis functions over the interval on which they are defined
  std::vector <double> bf_integrals_;
  // setup various stuff
  void setupBF();
  void setupInterval();
  void numericalBFIntegrals();
  virtual void setupDescription();
  virtual void setupBFIntegrals();
  template<typename T>
  void addKeywordToList(const std::string, const T);
  void addKeywordToList(const std::string, const bool);
public:
  static void registerKeywords(Keywords&);
  BasisFunctions(const ActionOptions&ao);
  bool hasBeenSet();
  std::string getType();
  std::string getDescription();
  unsigned int getOrder();
  unsigned int getNumberOfBasisFunctions();
  unsigned int getSize();
  bool arePeriodic();
  bool intervalBounded();
  double intervalMin();
  double intervalMax();
  double intervalRange();
  double intervalMean();
  double intervalDerivf();
  double getBasisFunctionIntegral(const unsigned int);
  std::vector<double> getBasisFunctionIntegrals();
  unsigned getNumberOfDerivatives();
  std::vector<std::string> getKeywordList();
  std::string getBasisFunctionDescription(const unsigned int);
  std::vector<std::string> getBasisFunctionDescriptions();
  //
  BasisFunctions();
  double translateArgument(const double, bool&);
  void apply();
  void calculate();
  // calculate the value for the n-th basis function
  virtual double getValue(const double, const unsigned int, double&, bool&)=0;
  // calcuate the values for all basis functions
  virtual void getAllValues(const double, double&, bool&, std::vector<double>&, std::vector<double>&)=0;
  //virtual void get2ndDerivaties(const double, std::vector<double>&)=0;
  void printInfo();
  std::string getKeywordString();
};

inline
bool BasisFunctions::hasBeenSet(){return has_been_set;}

inline
std::string BasisFunctions::getType(){return type_;}

inline
std::string BasisFunctions::getDescription(){return description_;}

inline
unsigned int BasisFunctions::getOrder(){return norder_;}

inline
unsigned int BasisFunctions::getNumberOfBasisFunctions(){return nbasis_;}

inline
unsigned int BasisFunctions::getSize(){return getNumberOfBasisFunctions();}

inline
bool BasisFunctions::arePeriodic(){return periodic_;}

inline
bool BasisFunctions::intervalBounded(){return interval_bounded_;}

inline
double BasisFunctions::intervalMin(){return interval_min_;}

inline
double BasisFunctions::intervalMax(){return interval_max_;}

inline
double BasisFunctions::intervalRange(){return interval_range_;}

inline
double BasisFunctions::intervalMean(){return interval_mean_;}

inline
double BasisFunctions::getBasisFunctionIntegral(unsigned int index){return bf_integrals_[index];}

inline
std::vector<double> BasisFunctions::getBasisFunctionIntegrals(){return bf_integrals_;}

inline
unsigned BasisFunctions::getNumberOfDerivatives(){return 0;}

inline
std::vector<std::string> BasisFunctions::getKeywordList(){return bf_keywords_;}

inline
std::string BasisFunctions::getBasisFunctionDescription(unsigned int index){return bf_description_[index];}

inline
std::vector<std::string> BasisFunctions::getBasisFunctionDescriptions(){return bf_description_;}


}

#endif

