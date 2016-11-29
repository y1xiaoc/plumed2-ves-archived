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
#include "tools/Tools.h"


namespace PLMD {

//+PLUMEDOC INTERNAL GAUSSIAN
/*
  Gaussian target distribution
*/
//+ENDPLUMEDOC

class TD_VonMises: public TargetDistribution {
  // properties of the Gaussians
  std::vector< std::vector<double> > sigmas_;
  std::vector< std::vector<double> > kappas_;
  std::vector< std::vector<double> > centers_;
  std::vector<double> weights_;
  std::vector<double> periods_;
  unsigned int ncenters_;
  double VonMisesDiagonal(const std::vector<double>&, const std::vector<double>&, const std::vector<double>&) const;
public:
  static void registerKeywords(Keywords&);
  explicit TD_VonMises(const TargetDistributionOptions& to);
  double getValue(const std::vector<double>&) const;
};


VES_REGISTER_TARGET_DISTRIBUTION(TD_VonMises,"VON_MISES")


void TD_VonMises::registerKeywords(Keywords& keys){
  TargetDistribution::registerKeywords(keys);
  keys.add("numbered","CENTER","The centers of the Von Mises distributions.");
  keys.add("numbered","SIGMA","The sigmas of the Von Mises distributions.");
  keys.add("optional","WEIGHTS","The weights of the Von Mises distributions.");
  keys.add("optional","PERIODS","The periods for each of the dimensions. By default they are 2*pi for each dimension.");
}


TD_VonMises::TD_VonMises( const TargetDistributionOptions& to ):
TargetDistribution(to),
sigmas_(0),
centers_(0),
weights_(0),
ncenters_(0)
{
  for(unsigned int i=1;; i++) {
    std::vector<double> tmp_center;
    if(!parseNumberedVector("CENTER",i,tmp_center) ){break;}
    centers_.push_back(tmp_center);
  }
  for(unsigned int i=1;; i++) {
    std::vector<double> tmp_sigma;
    if(!parseNumberedVector("SIGMA",i,tmp_sigma) ){break;}
    sigmas_.push_back(tmp_sigma);
  }
  if(centers_.size()==0 && sigmas_.size()==0){
    std::vector<double> tmp_center;
    if(parseVector("CENTER",tmp_center,true)){
      centers_.push_back(tmp_center);
    }
    std::vector<double> tmp_sigma;
    if(parseVector("SIGMA",tmp_sigma,true)){
      sigmas_.push_back(tmp_sigma);
    }
  }
  plumed_massert(centers_.size()==sigmas_.size(),"there has to be an equal amount of CENTER and SIGMA keywords");
  if(centers_.size()==0){
    plumed_merror("VON_MISES: CENTER and SIGMA keywords seem to be missing. Note that numbered keywords start at CENTER1 and SIGMA1.");
  }
  //
  setDimension(centers_[0].size());
  ncenters_ = centers_.size();
  //
  if(ncenters_>1){
    plumed_merror("For now VON_MISES only supports one center. Use LINEAR_COMBINATION to create a sum of von Mises distributions.");
  }
  // check centers and sigmas
  for(unsigned int i=0; i<ncenters_; i++) {
    plumed_massert(centers_[i].size()==getDimension(),"one of the CENTER keyword does not match the given dimension");
    plumed_massert(sigmas_[i].size()==getDimension(),"one of the CENTER keyword does not match the given dimension");
  }
  //
  kappas_.resize(sigmas_.size());
  for(unsigned int i=0; i<sigmas_.size();i++){
    kappas_[i].resize(sigmas_[i].size());
    for(unsigned int k=0; k<kappas_[i].size(); k++){
      kappas_[i][k] = 1.0/(sigmas_[i][k]*sigmas_[i][k]);
    }
  }
  //
  if(!parseVector("WEIGHTS",weights_,true)){weights_.assign(centers_.size(),1.0);}
  plumed_massert(centers_.size()==weights_.size(),"there has to be as many weights given in WEIGHTS as numbered CENTER keywords");
  //
  if(!parseVector("PERIODS",periods_,true)){periods_.assign(getDimension(),2*pi);}
  plumed_massert(periods_.size()==getDimension(),"the number of values given in PERIODS does not match the dimension of the distribution");
  //
  double sum_weights=0.0;
  for(unsigned int i=0;i<weights_.size();i++){sum_weights+=weights_[i];}
  for(unsigned int i=0;i<weights_.size();i++){weights_[i]/=sum_weights;}
  //
  setForcedNormalization();
  checkRead();
}


double TD_VonMises::getValue(const std::vector<double>& argument) const {
  double value=0.0;
  for(unsigned int i=0;i<ncenters_;i++){
    value+=weights_[i]*VonMisesDiagonal(argument, centers_[i], kappas_[i]);
  }
  return value;
}


double TD_VonMises::VonMisesDiagonal(const std::vector<double>& argument, const std::vector<double>& center, const std::vector<double>& kappa) const {
  double value = 1.0;
  for(unsigned int k=0; k<argument.size(); k++){
    double arg = kappa[k]*cos( ((2*pi)/periods_[k])*(argument[k]-center[k]) );
    value*=exp(arg);
  }
  return value;
}


}
