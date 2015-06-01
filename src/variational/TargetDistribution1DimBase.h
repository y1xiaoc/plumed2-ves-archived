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

class Grid;
class Keywords;

class TargetDistribution1DimOptions{
friend class TargetDistribution1DimRegister;
friend class TargetDistribution1DimBase;
private:
  std::vector<std::string> words;
public:
  TargetDistribution1DimOptions( const std::vector<std::string>& input);
};

class TargetDistribution1DimBase {
private:
/// Name of the one dimensional target distribution 
  std::string type;
/// The input to the target distribution
  std::vector<std::string> input;
/// is the target distribution normalize or not
  bool normalized;
protected:
/// Read a keywords from the input 
  template <class T>
  bool parse(const std::string& ,T& , bool optional=false);
/// Read a keywords vector from the input 
  template <class T>
  bool parseVector(const std::string& ,std::vector<T>& , bool optional=false);
/// Read a flag from the input
  void parseFlag(const std::string& key, bool& t);
public:
/// keywords
  static void registerKeywords( Keywords&);
  TargetDistribution1DimBase( const TargetDistribution1DimOptions& to );
  virtual ~TargetDistribution1DimBase();
/// Check everything was read in
  void checkRead() const ;
/// Return a description 
  std::string description();
/// Overwrite this to have a more descriptive output
  virtual std::string rest_of_description(){ return ""; };
/// is the target distribution normalize or not
   bool isNormalized()const{return normalized;};
/// set the that target distribution is normalized
   void setNormalized(){normalized=true;};
   void setNotNormalized(){normalized=false;};  
/// get type of distribution
   std::string getType()const{return type;};
/// calculate the target distribution itself
  virtual double distribution(const double)=0;
/// write the distribution out to file
  static void writeDistributionToFile(const std::string, const std::string, const double, const double, const unsigned int );
  void calculateDistributionOnGrid(Grid*);
};

template <class T>
bool TargetDistribution1DimBase::parse( const std::string& key, T& t, bool optional){
  bool found=Tools::parse(input,key,t);
  if(!optional && !found) plumed_merror("target distribution " + type + " requires " + key + " keyword");
  return found;
}

template <class T>
bool TargetDistribution1DimBase::parseVector( const std::string& key, std::vector<T>& t , bool optional){
  bool found=Tools::parseVector(input,key,t);
  if(!optional && !found) plumed_merror("target distribution " + type + " requires " + key + " keyword");
  return found;
}

}
#endif
