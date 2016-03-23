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
#include "TargetDistribution.h"
#include "TargetDistributionRegister.h"
#include "ves_biases/VesBias.h"

#include "core/Value.h"
#include "tools/Grid.h"
#include "tools/File.h"
#include "tools/Keywords.h"

namespace PLMD {

class MarginalWeight:public WeightBase{
  public:
    explicit MarginalWeight(){}
    double projectInnerLoop(double &input, double &v){return  input+v;}
    double projectOuterLoop(double &v){return v;}
};

TargetDistributionOptions::TargetDistributionOptions( const std::vector<std::string>& input):
words(input)
{
}


void TargetDistribution::registerKeywords( Keywords& keys ){
}


TargetDistribution::TargetDistribution( const TargetDistributionOptions& to):
type_(to.words[0]),
input(to.words),
normalized_(false),
dimension_(1),
action_pntr_(NULL),
vesbias_pntr_(NULL)
{
  input.erase( input.begin() );
}


TargetDistribution::~TargetDistribution() {
}


void TargetDistribution::setDimension(const unsigned int dimension){
  dimension_=dimension;
}


void TargetDistribution::linkVesBias(bias::VesBias* vesbias_pntr_in){
  vesbias_pntr_ = vesbias_pntr_in;
  action_pntr_ = static_cast<Action*>(vesbias_pntr_in);
}


void TargetDistribution::linkAction(Action* action_pntr_in){
  action_pntr_ = action_pntr_in;
}


void TargetDistribution::parseFlag(const std::string& key, bool& t) {
  Tools::parseFlag(input,key,t);
}


void TargetDistribution::checkRead() const {
  if(!input.empty()){
     std::string msg="cannot understand the following words from the target distribution input : ";
     for(unsigned i=0;i<input.size();++i) msg = msg + input[i] + ", ";
     plumed_merror(msg);
  }
}


std::string TargetDistribution::description() {
  std::string str="Type: " + type_;
  return str;
}


void TargetDistribution::writeDistributionToFile(const std::string& filepath, const std::vector<std::string>& min, const std::vector<std::string>& max, const std::vector<unsigned int>& nbins) {
  unsigned int dimension = getDimension();
  plumed_assert(min.size()==dimension);
  plumed_assert(max.size()==dimension);
  plumed_assert(nbins.size()==dimension);
  std::vector<Value*> arguments(dimension);
  for (unsigned int i=0; i < dimension; i++) {
    std::string is; Tools::convert(i+1,is);
    arguments[i]= new Value(NULL,"arg"+is,false);
    arguments[i]->setNotPeriodic();
  }
  Grid* grid_pntr = new Grid(getType(),arguments,min,max,nbins,false,false);
  //
  calculateDistributionOnGrid(grid_pntr);
  TargetDistribution::writeProbGridToFile(filepath,grid_pntr);

  // delete stuff
  delete grid_pntr;
  for (unsigned int i=0; i < dimension; i++) {delete arguments[i];}
}


void TargetDistribution::writeDistributionToFile(const std::string& filepath, const std::string& keywords, const std::vector<std::string>& min, const std::vector<std::string>& max, const std::vector<unsigned int>& nbins) {
  // create distribtion
  std::vector<std::string> words = Tools::getWords(keywords);
  TargetDistribution* targetdist_pntr=targetDistributionRegister().create( (words) );
  targetdist_pntr->writeDistributionToFile(filepath,min,max,nbins);
  delete targetdist_pntr;
}


void TargetDistribution::calculateDistributionOnGrid(Grid* grid_pntr) const {
  plumed_massert(grid_pntr->getDimension()==dimension_,"Grid is of the wrong dimension");
  for(unsigned int l=0; l<grid_pntr->getSize(); l++)
  {
   std::vector<double> argument=grid_pntr->getPoint(l);
   double value=getValue(argument);
   grid_pntr->setValue(l,value);
  }
}


void TargetDistribution::writeProbGridToFile(const std::string& filepath, Grid* grid_pntr, const bool do_projections) {
  OFile file;
  file.setBackupString("bck");
  file.open(filepath);
  grid_pntr->writeToFile(file);
  file.close();
  if(do_projections && grid_pntr->getDimension()>1){
    MarginalWeight* Pw = new MarginalWeight();
    std::vector<std::string> argnames = grid_pntr->getArgNames();
    std::vector<double> binDeltas = grid_pntr->getDx();
    for(unsigned int i=0; i<argnames.size(); i++){
      std::vector<std::string> arg(1,argnames[i]);
      Grid proj_grid = grid_pntr->project(arg,Pw);
      // scale with the bin volume used for the integral such that the
      // marginals are proberly normalized
      double intVol = 1.0;
      for(unsigned int l=0; l<binDeltas.size(); l++){
        if(l!=i){intVol*=binDeltas[l];}
      }
      proj_grid.scaleAllValuesAndDerivatives(intVol);
      //
      OFile file2;
      file2.setBackupString("bck");
      std::string proj_fname = argnames[i];
      proj_fname = FileBase::appendSuffix(filepath,".proj-"+proj_fname);
      file2.open(proj_fname);
      proj_grid.writeToFile(file2);
      file2.close();
    }
    delete Pw;
  }
}


void TargetDistribution::calculateSeperableDistributionOnGrid(Grid* grid_pntr, std::vector<TargetDistribution*> targetdist_pntrs) {
  unsigned int ntargetdist = targetdist_pntrs.size();
  plumed_massert(grid_pntr->getDimension()==ntargetdist,"dimension of Grid doesn't match the number of seperable one-dimensional target distribtions");
  for(unsigned int k=0; k<ntargetdist; k++){
    plumed_massert(targetdist_pntrs[k]!=NULL,"all the target distribtions must be defined");
    plumed_massert(targetdist_pntrs[k]->getDimension()==1,"all the target distribtions must be one-dimensional");
  }
  for(unsigned int l=0; l<grid_pntr->getSize(); l++)
  {
   std::vector<double> argument=grid_pntr->getPoint(l);
   double value=1;
   std::vector<double> arg1d(1);
   for(unsigned int k=0; k<ntargetdist; k++){
     arg1d[0] = argument[k];
     value*=targetdist_pntrs[k]->getValue(arg1d);
   }
   grid_pntr->setValue(l,value);
  }
}


}
