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

class ProductCombinationOfDistributions: public TargetDistribution {
  std::vector<TargetDistribution*> distributions_;
  unsigned int ndist_;
public:
  static void registerKeywords(Keywords&);
  explicit ProductCombinationOfDistributions(const TargetDistributionOptions& to);
  ~ProductCombinationOfDistributions();
  double getValue(const std::vector<double>&) const;
};


VES_REGISTER_TARGET_DISTRIBUTION(ProductCombinationOfDistributions,"PRODUCT_COMBINATION")


void ProductCombinationOfDistributions::registerKeywords(Keywords& keys){
  TargetDistribution::registerKeywords(keys);
  keys.add("numbered","ARG","The target distributions to be used in the product combination for each argument");
}


ProductCombinationOfDistributions::ProductCombinationOfDistributions( const TargetDistributionOptions& to ):
TargetDistribution(to),
distributions_(0),
ndist_(0)
{
  for(unsigned int i=1;; i++) {
    std::string keywords;
    if(!parseNumbered("ARG",i,keywords) ){break;}
    std::vector<std::string> words = Tools::getWords(keywords);
    TargetDistribution* dist_tmp = targetDistributionRegister().create( (words) );
    if(dist_tmp->getDimension()!=1){
      plumed_merror("PRODUCT_COMBINATION: all the target distributions should be one dimensional");
    }
    distributions_.push_back(dist_tmp);
  }
  ndist_ = distributions_.size();
  setDimension(ndist_);
  setNormalized();
  checkRead();
}


ProductCombinationOfDistributions::~ProductCombinationOfDistributions(){
  for(unsigned int i=0; i<ndist_; i++){
    delete distributions_[i];
  }
}


double ProductCombinationOfDistributions::getValue(const std::vector<double>& argument) const {
  double value=1.0;
  for(unsigned int i=0; i<ndist_; i++){
    std::vector<double> arg_tmp(1);
    arg_tmp[0]=argument[i];
    value *= distributions_[i]->getValue(arg_tmp);
  }
  return value;
}



}
