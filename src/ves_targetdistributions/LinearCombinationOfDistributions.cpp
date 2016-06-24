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

namespace PLMD {

//+PLUMEDOC INTERNAL GAUSSIAN
/*
  Gaussian target distribution
*/
//+ENDPLUMEDOC

class LinearCombinationOfDistributions: public TargetDistribution {
  std::vector<TargetDistribution*> distributions_;
  std::vector<double> weights_;
  unsigned int ndist_;
public:
  static void registerKeywords(Keywords&);
  explicit LinearCombinationOfDistributions(const TargetDistributionOptions& to);
  ~LinearCombinationOfDistributions();
  double getValue(const std::vector<double>&) const;
};


VES_REGISTER_TARGET_DISTRIBUTION(LinearCombinationOfDistributions,"LINEAR_COMBINATION")


void LinearCombinationOfDistributions::registerKeywords(Keywords& keys){
  TargetDistribution::registerKeywords(keys);
  keys.add("numbered","DISTRIBUTION","The target distributions to be used in the linear combination.");
  keys.add("optional","WEIGHTS","The weights of target distributions. If no weights are given the distributions are weighted equally.");
  keys.addFlag("DO_NOT_NORMALIZE",false,"If the weight should not be normalized.");
}


LinearCombinationOfDistributions::LinearCombinationOfDistributions( const TargetDistributionOptions& to ):
TargetDistribution(to),
distributions_(0),
weights_(0),
ndist_(0)
{
  for(unsigned int i=1;; i++) {
    std::string keywords;
    if(!parseNumbered("DISTRIBUTION",i,keywords) ){break;}
    std::vector<std::string> words = Tools::getWords(keywords);
    TargetDistribution* dist_tmp = targetDistributionRegister().create( (words) );
    distributions_.push_back(dist_tmp);
  }
  setDimension(distributions_[0]->getDimension());
  ndist_ = distributions_.size();

  for(unsigned int i=0; i<ndist_; i++){
    plumed_massert(getDimension()==distributions_[0]->getDimension(),"all distributions have to have the same dimension");
  }
  //
  if(!parseVector("WEIGHTS",weights_,true)){weights_.assign(distributions_.size(),1.0);}
  plumed_massert(distributions_.size()==weights_.size(),"there has to be as many weights given in WEIGHTS as numbered DISTRIBUTION keywords");
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
    delete distributions_[i];
  }
}


double LinearCombinationOfDistributions::getValue(const std::vector<double>& argument) const {
  double value=0.0;
  for(unsigned int i=0; i<ndist_; i++){
    value += weights_[i] * distributions_[i]->getValue(argument);
  }
  return value;
}



}
