/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2012-2014 The plumed team
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
#include "tools/Keywords.h"

namespace PLMD {

//+PLUMEDOC INTERNAL GAUSSIAN
/*
  Gaussian target distribution
*/
//+ENDPLUMEDOC

class LinearCombinationOfDistributions: public TargetDistributionBase {
  // properties of the Gaussians
  std::vector<TargetDistributionBase*> distributions;
  std::vector<double> weights;
  unsigned int ndist;
public:
  static void registerKeywords(Keywords&);
  LinearCombinationOfDistributions(const TargetDistributionOptions& to);
  ~LinearCombinationOfDistributions();
  double getValue(const std::vector<double>) const;
};


VARIATIONAL_REGISTER_TARGET_DISTRIBUTION(LinearCombinationOfDistributions,"LINEAR_COMBINATION")


void LinearCombinationOfDistributions::registerKeywords(Keywords& keys){
  TargetDistributionBase::registerKeywords(keys);
  keys.add("numbered","DISTRIBUTION","t.");
  keys.add("optional","WEIGHTS","The weights of the Gaussians.");
  keys.addFlag("DO_NOT_NORMALIZE",false,"If the weight should not be normalized.");
}


LinearCombinationOfDistributions::LinearCombinationOfDistributions( const TargetDistributionOptions& to ):
TargetDistributionBase(to)
{
  for(unsigned int i=0;; i++) {
    std::string keywords;
    if(!parseNumbered("DISTRIBUTION",i,keywords) ){break;}
    std::vector<std::string> words = Tools::getWords(keywords);
    TargetDistributionBase* dist_tmp = targetDistributionRegister().create( (words) );
    distributions.push_back(dist_tmp);
  }
  setDimension(distributions[0]->getDimension());
  ndist = distributions.size();

  for(unsigned int i=0; i<ndist; i++){
    plumed_massert(getDimension()==distributions[0]->getDimension(),"all distributions have to have the same dimension");
  }
  //
  if(!parseVector("WEIGHTS",weights,true)){weights.assign(distributions.size(),1.0);}
  plumed_massert(distributions.size()==weights.size(),"there has to be as many weights given in WEIGHTS as numbered DISTRIBUTION keywords");
  //
  bool do_not_normalize=false;
  parseFlag("DO_NOT_NORMALIZE",do_not_normalize);
  if(!do_not_normalize){
    double sum_weights=0.0;
    for(unsigned int i=0;i<weights.size();i++){sum_weights+=weights[i];}
    for(unsigned int i=0;i<weights.size();i++){weights[i]/=sum_weights;}
    setNormalized();
  }
  else{
    setNotNormalized();
  }
  checkRead();
}


LinearCombinationOfDistributions::~LinearCombinationOfDistributions(){
  for(unsigned int i=0; i<ndist; i++){
    delete distributions[i];
  }
}


double LinearCombinationOfDistributions::getValue(const std::vector<double> argument) const {
  double value=0.0;
  for(unsigned int i=0; i<ndist; i++){
    value += weights[i] * distributions[i]->getValue(argument);
  }
  return value;
}



}
