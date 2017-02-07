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


namespace PLMD{
namespace ves{

//+PLUMEDOC VES_TARGETDIST_HIDDEN EXPONENTIAL_POWER
/*
Target distribution given by a sum of exponential power distributions (static).

\par Examples

*/
//+ENDPLUMEDOC

class TD_ExponentialPower: public TargetDistribution {
  std::vector< std::vector<double> > centers_;
  std::vector< std::vector<double> > alphas_;
  std::vector< std::vector<double> > betas_;
  std::vector< std::vector<double> > normalization_;
  std::vector<double> weights_;
  unsigned int ncenters_;
  double ExponentialPowerDiagonal(const std::vector<double>&, const std::vector<double>&, const std::vector<double>&, const std::vector<double>&, const std::vector<double>&) const;
public:
  static void registerKeywords(Keywords&);
  explicit TD_ExponentialPower(const TargetDistributionOptions& to);
  double getValue(const std::vector<double>&) const;
};


VES_REGISTER_TARGET_DISTRIBUTION(TD_ExponentialPower,"EXPONENTIAL_POWER")


void TD_ExponentialPower::registerKeywords(Keywords& keys){
  TargetDistribution::registerKeywords(keys);
  keys.add("numbered","CENTER","The center of each exponential power distribution.");
  keys.add("numbered","ALPHA","The alpha parameters for each exponential power distribution.");
  keys.add("numbered","BETA","The beta parameters for each exponential power distribution.");
  keys.add("optional","WEIGHTS","The weights of the exponential power distribution. By default all are weighted equally.");
}


TD_ExponentialPower::TD_ExponentialPower( const TargetDistributionOptions& to ):
TargetDistribution(to),
centers_(0),
alphas_(0),
betas_(0),
normalization_(0),
weights_(0),
ncenters_(0)
{
  for(unsigned int i=1;; i++) {
    std::vector<double> tmp_center;
    if(!parseNumberedVector("CENTER",i,tmp_center) ){break;}
    centers_.push_back(tmp_center);
  }
  for(unsigned int i=1;; i++) {
    std::vector<double> tmp_alpha;
    if(!parseNumberedVector("ALPHA",i,tmp_alpha) ){break;}
    for(unsigned int k=0; k<tmp_alpha.size(); k++){
      if(tmp_alpha[k]<=0.0){plumed_merror(getName()+": the values given in ALPHA should be postive");}
    }
    alphas_.push_back(tmp_alpha);
  }
  for(unsigned int i=1;; i++) {
    std::vector<double> tmp_beta;
    if(!parseNumberedVector("BETA",i,tmp_beta) ){break;}
    for(unsigned int k=0; k<tmp_beta.size(); k++){
      if(tmp_beta[k]<=0.0){plumed_merror(getName()+": the values given in BETA should be postive");}
    }
    betas_.push_back(tmp_beta);
  }
  if(centers_.size()==0 && alphas_.size()==0 && betas_.size()==0){
    std::vector<double> tmp_center;
    if(parseVector("CENTER",tmp_center,true)){
      centers_.push_back(tmp_center);
    }
    std::vector<double> tmp_alpha;
    if(parseVector("ALPHA",tmp_alpha,true)){
      alphas_.push_back(tmp_alpha);
    }
    std::vector<double> tmp_beta;
    if(parseVector("BETA",tmp_beta,true)){
      betas_.push_back(tmp_beta);
    }
  }
  //
  if(centers_.size()==0){
    plumed_merror(getName()+": CENTER keywords seem to be missing. Note that numbered keywords start at CENTER1.");
  }
  //
  if(centers_.size()!=alphas_.size() || centers_.size()!=betas_.size() ){
    plumed_merror(getName()+": there has to be an equal amount of CENTER, ALPHA, and BETA keywords");
  }
  //
  setDimension(centers_[0].size());
  ncenters_ = centers_.size();
  //
  // check centers and sigmas
  for(unsigned int i=0; i<ncenters_; i++) {
    if(centers_[i].size()!=getDimension()){
      plumed_merror(getName()+": one of the CENTER keyword does not match the given dimension");
    }
    if(alphas_[i].size()!=getDimension()){
      plumed_merror(getName()+": one of the ALPHA keyword does not match the given dimension");
    }
    if(betas_[i].size()!=getDimension()){
      plumed_merror(getName()+": one of the BETA keyword does not match the given dimension");
    }
  }
  //
  if(!parseVector("WEIGHTS",weights_,true)){weights_.assign(centers_.size(),1.0);}
  if(centers_.size()!=weights_.size()){
    plumed_merror(getName()+": there has to be as many weights given in WEIGHTS as numbered CENTER keywords");
  }
  //
  double sum_weights=0.0;
  for(unsigned int i=0;i<weights_.size();i++){sum_weights+=weights_[i];}
  for(unsigned int i=0;i<weights_.size();i++){weights_[i]/=sum_weights;}
  //
  normalization_.resize(ncenters_);
  for(unsigned int i=0; i<ncenters_; i++) {
    normalization_[i].resize(getDimension());
    for(unsigned int k=0; k<getDimension(); k++){
      normalization_[i][k] = 0.5*betas_[i][k]/(alphas_[i][k]*tgamma(1.0/betas_[i][k]));
    }
  }
  checkRead();
}


double TD_ExponentialPower::getValue(const std::vector<double>& argument) const {
  double value=0.0;
  for(unsigned int i=0;i<ncenters_;i++){
    value+=weights_[i]*ExponentialPowerDiagonal(argument,centers_[i],alphas_[i],betas_[i],normalization_[i]);
  }
  return value;
}


double TD_ExponentialPower::ExponentialPowerDiagonal(const std::vector<double>& argument, const std::vector<double>& center, const std::vector<double>& alpha, const std::vector<double>& beta, const std::vector<double>& normalization) const {
  double value = 1.0;
  for(unsigned int k=0; k<argument.size(); k++){
    double arg=(std::abs(argument[k]-center[k]))/alpha[k];
    arg = pow(arg,beta[k]);
    value*=normalization[k]*exp(-arg);
  }
  return value;
}



}
}
