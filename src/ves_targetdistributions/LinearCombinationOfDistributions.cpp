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

#include "tools/Keywords.h"
#include "tools/Grid.h"

namespace PLMD {

//+PLUMEDOC INTERNAL GAUSSIAN
/*
  Gaussian target distribution
*/
//+ENDPLUMEDOC

class Grid;
class Action;

namespace bias{
  class VesBias;
}

class LinearCombinationOfDistributions: public TargetDistribution {
private:
  std::vector<TargetDistribution*> distribution_pntrs_;
  std::vector<Grid*> grid_pntrs_;
  std::vector<double> weights_;
  unsigned int ndist_;
  void setupAdditionalGrids(const std::vector<Value*>&, const std::vector<std::string>&, const std::vector<std::string>&, const std::vector<unsigned int>&);
public:
  static void registerKeywords(Keywords&);
  explicit LinearCombinationOfDistributions(const TargetDistributionOptions& to);
  void updateGrid();
  double getValue(const std::vector<double>&) const;
  ~LinearCombinationOfDistributions();
  //
  void linkVesBias(bias::VesBias*);
  void linkAction(Action*);
  void linkFesGrid(Grid*);
  void linkBiasGrid(Grid*);
  //
};


VES_REGISTER_TARGET_DISTRIBUTION(LinearCombinationOfDistributions,"LINEAR_COMBINATION")


void LinearCombinationOfDistributions::registerKeywords(Keywords& keys){
  TargetDistribution::registerKeywords(keys);
  keys.add("numbered","DISTRIBUTION","The target distributions to be used in the linear combination.");
  keys.add("optional","WEIGHTS","The weights of target distributions. If no weights are given the distributions are weighted equally.");
  keys.addFlag("DO_NOT_NORMALIZE",false,"If the weight should not be normalized. Be warned that this will lead to non-normalized distributions and most likely stange results.");
}


LinearCombinationOfDistributions::LinearCombinationOfDistributions( const TargetDistributionOptions& to ):
TargetDistribution(to),
distribution_pntrs_(0),
grid_pntrs_(0),
weights_(0),
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
  if(!parseVector("WEIGHTS",weights_,true)){weights_.assign(distribution_pntrs_.size(),1.0);}
  plumed_massert(distribution_pntrs_.size()==weights_.size(),"there has to be as many weights given in WEIGHTS as numbered DISTRIBUTION keywords");
  //
  bool do_not_normalize=false;
  parseFlag("DO_NOT_NORMALIZE",do_not_normalize);
  if(!do_not_normalize){
    double sum_weights=0.0;
    for(unsigned int i=0;i<weights_.size();i++){sum_weights+=weights_[i];}
    for(unsigned int i=0;i<weights_.size();i++){weights_[i]/=sum_weights;}
    setNormalized();
  }
  else{
    setNotNormalized();
  }
  checkRead();
}


LinearCombinationOfDistributions::~LinearCombinationOfDistributions(){
  for(unsigned int i=0; i<ndist_; i++){
    delete distribution_pntrs_[i];
  }
}


double LinearCombinationOfDistributions::getValue(const std::vector<double>& argument) const {
  plumed_merror("getValue not implemented for LinearCombinationOfDistributions");
  return 0.0;
}


void LinearCombinationOfDistributions::setupAdditionalGrids(const std::vector<Value*>& arguments, const std::vector<std::string>& min, const std::vector<std::string>& max, const std::vector<unsigned int>& nbins) {
  for(unsigned int i=0; i<ndist_; i++){
    distribution_pntrs_[i]->setupGrids(arguments,min,max,nbins);
    if(distribution_pntrs_[i]->getDimension()!=this->getDimension()){
      plumed_merror("Error in LINEAR_COMBINATION: all target distribution need to have the same dimension");
    }
    grid_pntrs_[i]=distribution_pntrs_[i]->getTargetDistGridPntr();
  }
}


void LinearCombinationOfDistributions::updateGrid(){
  for(unsigned int i=0; i<ndist_; i++){
    distribution_pntrs_[i]->updateGrid();
  }
  for(Grid::index_t l=0; l<targetDistGrid().getSize(); l++){
    double value = 0.0;
    for(unsigned int i=0; i<ndist_; i++){
      value += weights_[i]*grid_pntrs_[i]->getValue(l);
    }
    targetDistGrid().setValue(l,value);
    logTargetDistGrid().setValue(l,-std::log(value));
  }
  logTargetDistGrid().setMinToZero();
}


void LinearCombinationOfDistributions::linkVesBias(bias::VesBias* vesbias_pntr_in){
  TargetDistribution::linkVesBias(vesbias_pntr_in);
  for(unsigned int i=0; i<ndist_; i++){
    distribution_pntrs_[i]->linkVesBias(vesbias_pntr_in);
  }
}


void LinearCombinationOfDistributions::linkAction(Action* action_pntr_in){
  TargetDistribution::linkAction(action_pntr_in);
  for(unsigned int i=0; i<ndist_; i++){
    distribution_pntrs_[i]->linkAction(action_pntr_in);
  }
}


void LinearCombinationOfDistributions::linkBiasGrid(Grid* bias_grid_pntr_in){
  TargetDistribution::linkBiasGrid(bias_grid_pntr_in);
  for(unsigned int i=0; i<ndist_; i++){
    distribution_pntrs_[i]->linkBiasGrid(bias_grid_pntr_in);
  }
}


void LinearCombinationOfDistributions::linkFesGrid(Grid* fes_grid_pntr_in){
  TargetDistribution::linkFesGrid(fes_grid_pntr_in);
  for(unsigned int i=0; i<ndist_; i++){
    distribution_pntrs_[i]->linkFesGrid(fes_grid_pntr_in);
  }
}


}
