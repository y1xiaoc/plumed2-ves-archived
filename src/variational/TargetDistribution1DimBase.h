/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2013-2015 The plumed team
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
#ifndef __PLUMED_variational_TargetDistribution1DimBase_h
#define __PLUMED_variational_TargetDistribution1DimBase_h

#include <vector>
#include <string>
#include "tools/Exception.h"
#include "tools/Tools.h"

namespace PLMD {
class BasisFunctions;
namespace variational {

class TargetDistribution1DimOptions{
friend class TargetDistribution1DimRegister;
friend class TargetDistribution1DimBase;
private:
  std::vector<std::string> words;
  BasisFunctions* basisf;
public:
  TargetDistribution1DimOptions( const std::vector<std::string>& input, BasisFunctions* mybasisf );
};

class TargetDistribution1DimBase {
private:
/// Name of the one dimensional target distribution 
  std::string type;
/// The input to the target distribution
  std::vector<std::string> input;
/// A pointer to the underlying BasisFunction
  BasisFunctions* basisf;
/// is the target distribution normalize or not
  bool normalized;
  double normalization_factor;
protected:
/// Read a keywords from the input 
  template <class T>
  void parse(const std::string& ,T& );
/// Read a flag from the input
  void parseFlag(const std::string& key, bool& t);
   void setNormalizationFactor(double);
   double getNormalizationFactor(){return normalization_factor;}
public:
  TargetDistribution1DimBase( const TargetDistribution1DimOptions& to );
  virtual ~TargetDistribution1DimBase();
/// Check everything was read in
  void checkRead() const ;
/// Return a description of the landmark selection protocol
  std::string description();
/// Overwrite this to have a more descriptive output
  virtual std::string rest_of_description(){ return ""; };
/// is the target distribution normalize or not
   bool isNormalized(){return normalized;};
   BasisFunctions* BF_pointer(){return basisf;};
/// calculate the target distribution itself
  virtual double calculate_ps(double)=0;
};

template <class T>
void TargetDistribution1DimBase::parse( const std::string& key, T& t ){
  bool found=Tools::parse(input,key,t);
  if(!found) plumed_merror("target distribution " + type + " requires " + key + " keyword");
}

}
}
#endif
