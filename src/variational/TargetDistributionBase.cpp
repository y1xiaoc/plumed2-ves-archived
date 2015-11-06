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
#include "TargetDistributionBase.h"
#include "TargetDistributionRegister.h"
#include "core/Value.h"
#include "tools/Grid.h"
#include "tools/File.h"
#include "tools/Keywords.h"

namespace PLMD {

TargetDistributionOptions::TargetDistributionOptions( const std::vector<std::string>& input):
words(input)
{
}


void TargetDistributionBase::registerKeywords( Keywords& keys ){
  keys.add("optional","DIMENSION","The dimension of the target distribtion, is by default 1.");
}


TargetDistributionBase::TargetDistributionBase( const TargetDistributionOptions& to):
type(to.words[0]),
input(to.words),
normalized_(false),
dimension_(1)
{
  input.erase( input.begin() );
  parse("DIMENSION",dimension_,true);
}


TargetDistributionBase::~TargetDistributionBase() {
}


void TargetDistributionBase::parseFlag(const std::string& key, bool& t) {
  Tools::parseFlag(input,key,t);
}


void TargetDistributionBase::checkRead() const {
  if(!input.empty()){
     std::string msg="cannot understand the following words from the target distribution input : ";
     for(unsigned i=0;i<input.size();++i) msg = msg + input[i] + ", ";
     plumed_merror(msg);
  }
}


std::string TargetDistributionBase::description() {
  std::string str="Type: " + type;
  return str;
}


void TargetDistributionBase::writeDistributionToFile(const std::string filepath, const std::string keywords, const std::vector<std::string> min, const std::vector<std::string> max, const std::vector<unsigned int>  nbins) {
  // create distribtion
  std::vector<std::string> words = Tools::getWords(keywords);
  TargetDistributionBase* distribution=targetDistributionRegister().create( (words) );
  unsigned int dimension = distribution->getDimension();
  // create grid
  std::vector<Value*> arguments(dimension);
  for (unsigned int i=0; i < dimension; i++) {
    arguments[i]= new Value(NULL,"argument",false);
    arguments[i]->setNotPeriodic();
  }
  Grid* grid = new Grid(distribution->getType(),arguments,min,max,nbins,false,false);
  //
  distribution->calculateDistributionOnGrid(grid);
  // write to file
  OFile file; file.open(filepath);
  grid->writeToFile(file);
  file.close();
  // delete stuff
  delete grid;
  delete distribution;
  for (unsigned int i=0; i < dimension; i++) {delete arguments[i];}
}


void TargetDistributionBase::calculateDistributionOnGrid(Grid* grid_ptr){
  plumed_massert(grid_ptr->getDimension()==dimension_,"Grid is if the wrong dimension");
  for(unsigned int l=0; l<grid_ptr->getSize(); l++)
  {
   std::vector<double> argument=grid_ptr->getPoint(l);
   double value=distribution(argument);
   grid_ptr->setValue(l,value);
  }
}


}
