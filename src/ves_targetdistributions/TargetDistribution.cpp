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
#include "TargetDistribution.h"
#include "TargetDistributionRegister.h"
#include "ves_biases/VesBias.h"

#include "core/Value.h"
#include "tools/Grid.h"
#include "ves_tools/GridProjWeights.h"
#include "ves_tools/GridIntegrationWeights.h"
#include "tools/File.h"
#include "tools/Keywords.h"

namespace PLMD {

TargetDistributionOptions::TargetDistributionOptions( const std::vector<std::string>& input):
words(input){}


void TargetDistribution::registerKeywords( Keywords& keys ){
}


TargetDistribution::TargetDistribution( const TargetDistributionOptions& to):
name_(to.words[0]),
input(to.words),
type_(static_targetdist),
normalized_(false),
dimension_(0),
targetdist_grid_pntr_(NULL),
log_targetdist_grid_pntr_(NULL),
action_pntr_(NULL),
vesbias_pntr_(NULL),
needs_bias_grid_(false),
needs_fes_grid_(false),
fes_grid_pntr_(NULL),
bias_grid_pntr_(NULL),
static_grid_calculated(false)
{
  input.erase( input.begin() );
}


TargetDistribution::~TargetDistribution() {
  if(targetdist_grid_pntr_!=NULL){
    delete targetdist_grid_pntr_;
  }
  if(log_targetdist_grid_pntr_!=NULL){
    delete log_targetdist_grid_pntr_;
  }
}


double TargetDistribution::getBeta() const {
  plumed_massert(vesbias_pntr_!=NULL,"The VesBias has to be linked to use TargetDistribution::getBeta");
  return vesbias_pntr_->getBeta();
}


void TargetDistribution::setDimension(const unsigned int dimension){
  plumed_massert(dimension_==0,"setDimension: the dimension of the target distribution has already been set!");
  dimension_=dimension;
}


void TargetDistribution::linkVesBias(bias::VesBias* vesbias_pntr_in){
  vesbias_pntr_ = vesbias_pntr_in;
  action_pntr_ = static_cast<Action*>(vesbias_pntr_in);
}


void TargetDistribution::linkAction(Action* action_pntr_in){
  action_pntr_ = action_pntr_in;
}


void TargetDistribution::linkBiasGrid(Grid* bias_grid_pntr_in){
  bias_grid_pntr_=bias_grid_pntr_in;
}


void TargetDistribution::linkFesGrid(Grid* fes_grid_pntr_in){
  fes_grid_pntr_=fes_grid_pntr_in;
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
  std::string str="Type: " + name_;
  return str;
}


void TargetDistribution::setupGrids(const std::vector<Value*>& arguments, const std::vector<std::string>& min, const std::vector<std::string>& max, const std::vector<unsigned int>& nbins) {
  if(getDimension()==0){
    setDimension(arguments.size());
  }
  unsigned int dimension = getDimension();
  plumed_massert(arguments.size()==dimension,"TargetDistribution::setupGrids: mismatch between number of values given for grid parameters");
  plumed_massert(min.size()==dimension,"TargetDistribution::setupGrids: mismatch between number of values given for grid parameters");
  plumed_massert(max.size()==dimension,"TargetDistribution::setupGrids: mismatch between number of values given for grid parameters");
  plumed_massert(nbins.size()==dimension,"TargetDistribution::setupGrids: mismatch between number of values given for grid parameters");
  targetdist_grid_pntr_ =     new Grid("targetdist",arguments,min,max,nbins,false,false);
  log_targetdist_grid_pntr_ = new Grid("log_targetdist",arguments,min,max,nbins,false,false);
  setupAdditionalGrids(arguments,min,max,nbins);
}


void TargetDistribution::calculateStaticDistributionGrid(){
  if(static_grid_calculated){return;}
  plumed_massert(isStatic(),"this should only be used for static distributions");
  plumed_massert(targetdist_grid_pntr_!=NULL,"the grids have not been setup using setupGrids!!");
  plumed_massert(log_targetdist_grid_pntr_!=NULL,"the grids have not been setup using setupGrids!!!!");
  for(Grid::index_t l=0; l<targetdist_grid_pntr_->getSize(); l++)
  {
   std::vector<double> argument = targetdist_grid_pntr_->getPoint(l);
   double value = getValue(argument);
   targetdist_grid_pntr_->setValue(l,value);
   log_targetdist_grid_pntr_->setValue(l,-std::log(value));
  }
  log_targetdist_grid_pntr_->setMinToZero();
  static_grid_calculated = true;
}


double TargetDistribution::integrateGrid(const Grid* grid_pntr){
  std::vector<double> integration_weights = GridIntegrationWeights::getIntegrationWeights(grid_pntr);
  double sum = 0.0;
  for(Grid::index_t l=0; l<grid_pntr->getSize(); l++){
    sum += integration_weights[l]*grid_pntr->getValue(l);
  }
  return sum;
}


double TargetDistribution::normalizeGrid(Grid* grid_pntr){
  double normalization = TargetDistribution::integrateGrid(grid_pntr);
  grid_pntr->scaleAllValuesAndDerivatives(1.0/normalization);
  return normalization;
}


Grid TargetDistribution::getMarginalDistributionGrid(Grid* grid_pntr, const std::vector<std::string>& args) {
  plumed_massert(grid_pntr->getDimension()>1,"doesn't make sense calculating the marginal distribution for a one-dimensional distribution");
  plumed_massert(args.size()<grid_pntr->getDimension(),"the number of arguments for the marginal distribution should be less than the dimension of the full distribution");
  //
  std::vector<std::string> argnames = grid_pntr->getArgNames();
  std::vector<unsigned int> args_index(0);
  for(unsigned int i=0; i<argnames.size(); i++){
    for(unsigned int l=0; l<args.size(); l++){
      if(argnames[i]==args[l]){args_index.push_back(i);}
    }
  }
  plumed_massert(args.size()==args_index.size(),"getMarginalDistributionGrid: problem with the arguments of the marginal");
  //
  MarginalWeight* Pw = new MarginalWeight();
  Grid proj_grid = grid_pntr->project(args,Pw);
  delete Pw;
  //
  // scale with the bin volume used for the integral such that the
  // marginals are proberly normalized to 1.0
  double intVol = grid_pntr->getBinVolume();
  for(unsigned int l=0; l<args_index.size(); l++){
    intVol/=grid_pntr->getDx()[args_index[l]];
  }
  proj_grid.scaleAllValuesAndDerivatives(intVol);
  //
  return proj_grid;
}


Grid TargetDistribution::getMarginal(const std::vector<std::string>& args){
  return TargetDistribution::getMarginalDistributionGrid(targetdist_grid_pntr_,args);
}


}
