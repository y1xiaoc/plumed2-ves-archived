/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2015-2016 The ves-code team
   (see the PEOPLE-VES file at the root of the distribution for a list of names)

   See http://www.ves-code.org for more information.

   This file is part of ves-code, version 1.

   ves-code is free software: you can redistribute it and/or modify
   it under the terms of the GNU Lesser General Public License as published by
   the Free Software Foundation, either version 3 of the License, or
   (at your option) any later version.

   ves-code is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU Lesser General Public License for more details.

   You should have received a copy of the GNU Lesser General Public License
   along with ves-code.  If not, see <http://www.gnu.org/licenses/>.
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
#ifndef __PLUMED_ves_targetdistributions_TargetDistribution_h
#define __PLUMED_ves_targetdistributions_TargetDistribution_h

#include <vector>
#include <string>
#include "tools/Exception.h"
#include "tools/Tools.h"

namespace PLMD {

class Grid;
class Value;
class Keywords;
class Action;

namespace bias{
  class VesBias;
}

class TargetDistributionOptions{
friend class TargetDistributionRegister;
friend class TargetDistribution;
private:
  std::vector<std::string> words;
public:
  TargetDistributionOptions( const std::vector<std::string>& input);
};

class TargetDistribution {
private:
  // Name of the target distribution
  std::string name_;
  // The input to the target distribution
  std::vector<std::string> input;
  enum TargetDistType {
    static_targetdist,
    dynamic_targetdist
  } type_;
  // is the target distribution normalized
  bool normalized_;
  // dimension of the distribution
  unsigned int dimension_;
  //
  Grid* targetdist_grid_pntr_;
  Grid* log_targetdist_grid_pntr_;
  //
  Action* action_pntr_;
  bias::VesBias* vesbias_pntr_;
  //
  void calculateStaticDistributionGrid();
protected:
  // Read a keywords from the input
  template <class T>
  bool parse(const std::string& ,T& , bool optional=false);
  template <class T>
  bool parseNumbered(const std::string& ,const unsigned int, T&);
  // Read a keywords vector from the input
  template <class T>
  bool parseVector(const std::string& ,std::vector<T>& , bool optional=false);
  template <class T>
  bool parseNumberedVector(const std::string& ,const unsigned int, std::vector<T>&);
  // Read a flag from the input
  void parseFlag(const std::string& key, bool& t);
  //
  void setStatic(){type_=static_targetdist;}
  void setDynamic(){type_=dynamic_targetdist;}
  // set the that target distribution is normalized
  void setNormalized(){normalized_=true;};
  void setNotNormalized(){normalized_=false;};
  //
  bias::VesBias* getPntrToVesBias() const;
  Action* getPntrToAction() const;
  //
  virtual void setupAdditionalGrids(const std::vector<Value*>&, const std::vector<std::string>&, const std::vector<std::string>&, const std::vector<unsigned int>&) {}
  //
  void normalizeTargetDistGrid();
  //
  Grid& targetDistGrid() const {return *targetdist_grid_pntr_;}
  Grid& logTargetDistGrid() const {return *log_targetdist_grid_pntr_;}
public:
  // keywords
  static void registerKeywords( Keywords&);
  explicit TargetDistribution( const TargetDistributionOptions& to );
  virtual ~TargetDistribution();
  // Check everything was read in
  void checkRead() const ;
  // Return a description
  std::string description();
  // Overwrite this to have a more descriptive output
  virtual std::string rest_of_description(){ return ""; };
  //
  bool isStatic() const {return type_==static_targetdist;}
  bool isDynamic() const {return type_==dynamic_targetdist;}
  // is the target distribution normalize or not
  bool isNormalized() const {return normalized_;};
  //
  void setDimension(const unsigned int dimension);
  unsigned getDimension() const {return dimension_;}
  // get type of distribution
  std::string getName()const{return name_;};
  //
  void linkVesBias(bias::VesBias*);
  void linkAction(Action*);
  //
  Grid* getTargetDistGridPntr() const {return targetdist_grid_pntr_;}
  Grid* getLogTargetDistGridPntr() const {return log_targetdist_grid_pntr_;}
  // calculate the target distribution itself
  virtual double getValue(const std::vector<double>&) const = 0;
  //
  void setupGrids(const std::vector<Value*>&, const std::vector<std::string>&, const std::vector<std::string>&, const std::vector<unsigned int>&);
  //
  Grid getMarginal(const std::vector<std::string>&);
  //
  virtual void updateGrid() {calculateStaticDistributionGrid();};
  //
  static double integrateGrid(const Grid*);
  static double normalizeGrid(Grid*);
  static Grid getMarginalDistributionGrid(Grid*, const std::vector<std::string>&);
};


inline
bias::VesBias* TargetDistribution::getPntrToVesBias() const {
  plumed_massert(vesbias_pntr_!=NULL,"the VES bias has not been linked");
  return vesbias_pntr_;
}


inline
Action* TargetDistribution::getPntrToAction() const {
  plumed_massert(action_pntr_!=NULL,"the action has not been linked");
  return action_pntr_;
}


inline
void TargetDistribution::normalizeTargetDistGrid(){
  normalizeGrid(targetdist_grid_pntr_);
}


template <class T>
bool TargetDistribution::parse( const std::string& key, T& t, bool optional){
  bool found=Tools::parse(input,key,t);
  if(!optional && !found) plumed_merror("target distribution " + name_ + " requires " + key + " keyword");
  return found;
}


template<class T>
bool TargetDistribution::parseNumbered(const std::string&key, const unsigned int no, T&t) {
  std::string num; Tools::convert(no,num);
  return Tools::parse(input,key+num,t);
}


template <class T>
bool TargetDistribution::parseVector( const std::string& key, std::vector<T>& t , bool optional){
  bool found=Tools::parseVector(input,key,t);
  if(!optional && !found) plumed_merror("target distribution " + name_ + " requires " + key + " keyword");
  return found;
}


template <class T>
bool TargetDistribution::parseNumberedVector( const std::string& key, const unsigned int no, std::vector<T>& t) {
  std::string num; Tools::convert(no,num);
  return Tools::parseVector(input,key+num,t);
}




}
#endif
