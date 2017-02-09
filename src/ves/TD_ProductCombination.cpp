/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2016-2017 The ves-code team
   (see the PEOPLE-VES file at the root of this folder for a list of names)

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

#include "tools/Keywords.h"
#include "tools/Grid.h"
#include "GridIntegrationWeights.h"


namespace PLMD{

// class Grid;
class Action;

namespace ves{

//+PLUMEDOC VES_TARGETDIST PRODUCT_COMBINATION
/*
Target distribution given by product combination of distributions (static or dynamic).

\par Examples

*/
//+ENDPLUMEDOC

class VesBias;

class TD_ProductCombination: public TargetDistribution {
private:
  std::vector<TargetDistribution*> distribution_pntrs_;
  std::vector<Grid*> grid_pntrs_;
  unsigned int ndist_;
  void setupAdditionalGrids(const std::vector<Value*>&, const std::vector<std::string>&, const std::vector<std::string>&, const std::vector<unsigned int>&);
public:
  static void registerKeywords(Keywords&);
  explicit TD_ProductCombination(const TargetDistributionOptions& to);
  void updateGrid();
  double getValue(const std::vector<double>&) const;
  ~TD_ProductCombination();
  //
  void linkVesBias(VesBias*);
  void linkAction(Action*);
  //
  void linkBiasGrid(Grid*);
  void linkBiasWithoutCutoffGrid(Grid*);
  void linkFesGrid(Grid*);
  //
};


VES_REGISTER_TARGET_DISTRIBUTION(TD_ProductCombination,"PRODUCT_COMBINATION")


void TD_ProductCombination::registerKeywords(Keywords& keys){
  TargetDistribution::registerKeywords(keys);
  keys.add("numbered","DISTRIBUTION","The target distributions to be used in the product combination.");
  keys.use("BIAS_CUTOFF");
  keys.use("WELLTEMPERED_FACTOR");
  keys.use("SHIFT_TO_ZERO");
  keys.use("NORMALIZE");
}


TD_ProductCombination::TD_ProductCombination( const TargetDistributionOptions& to ):
TargetDistribution(to),
distribution_pntrs_(0),
grid_pntrs_(0),
ndist_(0)
{
  for(unsigned int i=1;; i++) {
    std::string keywords;
    if(!parseNumbered("DISTRIBUTION",i,keywords) ){break;}
    std::vector<std::string> words = Tools::getWords(keywords);
    TargetDistribution* dist_pntr_tmp = targetDistributionRegister().create( (words) );
    if(dist_pntr_tmp->isDynamic()){setDynamic();}
    if(dist_pntr_tmp->fesGridNeeded()){setFesGridNeeded();}
    if(dist_pntr_tmp->biasGridNeeded()){setBiasGridNeeded();}
    distribution_pntrs_.push_back(dist_pntr_tmp);
  }
  ndist_ = distribution_pntrs_.size();
  grid_pntrs_.assign(ndist_,NULL);
  //
  checkRead();
}


TD_ProductCombination::~TD_ProductCombination(){
  for(unsigned int i=0; i<ndist_; i++){
    delete distribution_pntrs_[i];
  }
}


double TD_ProductCombination::getValue(const std::vector<double>& argument) const {
  plumed_merror("getValue not implemented for TD_ProductCombination");
  return 0.0;
}


void TD_ProductCombination::setupAdditionalGrids(const std::vector<Value*>& arguments, const std::vector<std::string>& min, const std::vector<std::string>& max, const std::vector<unsigned int>& nbins) {
  for(unsigned int i=0; i<ndist_; i++){
    distribution_pntrs_[i]->setupGrids(arguments,min,max,nbins);
    if(distribution_pntrs_[i]->getDimension()!=this->getDimension()){
      plumed_merror(getName() + ": all target distribution must have the same dimension");
    }
    grid_pntrs_[i]=distribution_pntrs_[i]->getTargetDistGridPntr();
  }
}


void TD_ProductCombination::updateGrid(){
  for(unsigned int i=0; i<ndist_; i++){
    distribution_pntrs_[i]->update();
  }
  std::vector<double> integration_weights = GridIntegrationWeights::getIntegrationWeights(getTargetDistGridPntr());
  double norm = 0.0;
  for(Grid::index_t l=0; l<targetDistGrid().getSize(); l++){
    double value = 1.0;
    for(unsigned int i=0; i<ndist_; i++){
      value *= grid_pntrs_[i]->getValue(l);
    }
    norm += integration_weights[l]*value;
    targetDistGrid().setValue(l,value);
    logTargetDistGrid().setValue(l,-std::log(value));
  }
  targetDistGrid().scaleAllValuesAndDerivatives(1.0/norm);
  logTargetDistGrid().setMinToZero();
}


void TD_ProductCombination::linkVesBias(VesBias* vesbias_pntr_in){
  TargetDistribution::linkVesBias(vesbias_pntr_in);
  for(unsigned int i=0; i<ndist_; i++){
    distribution_pntrs_[i]->linkVesBias(vesbias_pntr_in);
  }
}


void TD_ProductCombination::linkAction(Action* action_pntr_in){
  TargetDistribution::linkAction(action_pntr_in);
  for(unsigned int i=0; i<ndist_; i++){
    distribution_pntrs_[i]->linkAction(action_pntr_in);
  }
}


void TD_ProductCombination::linkBiasGrid(Grid* bias_grid_pntr_in){
  TargetDistribution::linkBiasGrid(bias_grid_pntr_in);
  for(unsigned int i=0; i<ndist_; i++){
    distribution_pntrs_[i]->linkBiasGrid(bias_grid_pntr_in);
  }
}


void TD_ProductCombination::linkBiasWithoutCutoffGrid(Grid* bias_withoutcutoff_grid_pntr_in){
  TargetDistribution::linkBiasWithoutCutoffGrid(bias_withoutcutoff_grid_pntr_in);
  for(unsigned int i=0; i<ndist_; i++){
    distribution_pntrs_[i]->linkBiasWithoutCutoffGrid(bias_withoutcutoff_grid_pntr_in);
  }
}


void TD_ProductCombination::linkFesGrid(Grid* fes_grid_pntr_in){
  TargetDistribution::linkFesGrid(fes_grid_pntr_in);
  for(unsigned int i=0; i<ndist_; i++){
    distribution_pntrs_[i]->linkFesGrid(fes_grid_pntr_in);
  }
}


}
}
