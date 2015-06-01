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
#include "TargetDistribution1DimBase.h"
#include "TargetDistribution1DimRegister.h"
#include "core/Value.h"
#include "tools/Grid.h"
#include "tools/File.h"
#include "tools/Keywords.h"

namespace PLMD {

TargetDistribution1DimOptions::TargetDistribution1DimOptions( const std::vector<std::string>& input):
words(input)
{
}

void TargetDistribution1DimBase::registerKeywords( Keywords& keys ){
}

TargetDistribution1DimBase::TargetDistribution1DimBase( const TargetDistribution1DimOptions& to ):
type(to.words[0]),
input(to.words),
normalized(false)
{
  input.erase( input.begin() );
}

TargetDistribution1DimBase::~TargetDistribution1DimBase(){
}


void TargetDistribution1DimBase::parseFlag(const std::string& key, bool& t){
  Tools::parseFlag(input,key,t);
}

void TargetDistribution1DimBase::checkRead() const {
  if(!input.empty()){
     std::string msg="cannot understand the following words from landmark selection input : ";
     for(unsigned i=0;i<input.size();++i) msg = msg + input[i] + ", ";
     plumed_merror(msg); 
  }
}

std::string TargetDistribution1DimBase::description(){
  std::string str="Type: " + type;
  return str;
}

void TargetDistribution1DimBase::writeDistributionToFile(const std::string filepath, const std::string keywords, const double interval_min, const double interval_max, const unsigned int number_of_bins){
  // create distribtion 
  std::vector<std::string> words = Tools::getWords(keywords);
  TargetDistribution1DimBase* distribution=targetDistribution1DimRegister().create( (words) );
  // create grid 
  Value* arg = new Value(NULL,"argument",false);
  arg->setNotPeriodic();
  std::vector<Value*> arguments(1); arguments[0]=arg;
  std::vector<std::string> min(1); Tools::convert(interval_min,min[0]);
  std::vector<std::string> max(1); Tools::convert(interval_max,max[0]);
  std::vector<unsigned int> nbins(1); nbins[0]=number_of_bins;
  Grid* grid = new Grid(distribution->getType(),arguments,min,max,nbins,false,false);
  // 
  distribution->calculateDistributionOnGrid(grid);
  // write to file
  OFile file; file.open(filepath);
  grid->writeToFile(file);
  file.close();
  // delete stuff
  delete arg;
  delete grid;
  delete distribution;
}

void TargetDistribution1DimBase::calculateDistributionOnGrid(Grid* grid_ptr){
  plumed_massert(grid_ptr->getDimension()==1,"Grid should be one-dimensional");
  for(unsigned int l=0; l<grid_ptr->getSize(); l++)
  {
   std::vector<double> argument=grid_ptr->getPoint(l);
   double value=distribution(argument[0]);
   grid_ptr->setValue(l,value);
  }
}  


}

