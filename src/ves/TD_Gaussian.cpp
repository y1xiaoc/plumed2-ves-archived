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


namespace PLMD{
namespace ves{

//+PLUMEDOC VES_TARGETDIST GAUSSIAN
/*
Target distribution given by a sum of Gaussians (static).

\par Examples

*/
//+ENDPLUMEDOC

class TD_Gaussian: public TargetDistribution {
  std::vector< std::vector<double> > centers_;
  std::vector< std::vector<double> > sigmas_;
  std::vector< std::vector<double> > correlation_;
  std::vector<double> weights_;
  bool diagonal_;
  unsigned int ncenters_;
  double GaussianDiagonal(const std::vector<double>&, const std::vector<double>&, const std::vector<double>&, const bool normalize=true) const;
  double Gaussian2D(const std::vector<double>&, const std::vector<double>&, const std::vector<double>&, const std::vector<double>&, const bool normalize=true) const;
public:
  static void registerKeywords(Keywords&);
  explicit TD_Gaussian(const TargetDistributionOptions& to);
  double getValue(const std::vector<double>&) const;
};


VES_REGISTER_TARGET_DISTRIBUTION(TD_Gaussian,"GAUSSIAN")


void TD_Gaussian::registerKeywords(Keywords& keys){
  TargetDistribution::registerKeywords(keys);
  keys.add("numbered","CENTER","The center of each Gaussian distribution.");
  keys.add("numbered","SIGMA","The sigma parameters for each Gaussian distribution.");
  keys.add("numbered","CORRELATION","The correlation between the arguments for each Gaussian distribution. Currently only works for two-dimensional Gaussians.");
  keys.add("optional","WEIGHTS","The weights of the Gaussian distribution. By default all are weighted equally.");
}


TD_Gaussian::TD_Gaussian( const TargetDistributionOptions& to ):
TargetDistribution(to),
centers_(0),
sigmas_(0),
correlation_(0),
weights_(0),
diagonal_(true),
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
  //
  if(centers_.size()==0){
    plumed_merror(getName()+": CENTER keywords seem to be missing. Note that numbered keywords start at CENTER1.");
  }
  //
  if(centers_.size()!=sigmas_.size()){
    plumed_merror(getName()+": there has to be an equal amount of CENTER and SIGMA keywords");
  }
  //
  setDimension(centers_[0].size());
  ncenters_ = centers_.size();
  // check centers and sigmas
  for(unsigned int i=0; i<ncenters_; i++) {
    if(centers_[i].size()!=getDimension()){
      plumed_merror(getName()+": one of the CENTER keyword does not match the given dimension");
    }
    if(sigmas_[i].size()!=getDimension()){
      plumed_merror(getName()+": one of the SIGMA keyword does not match the given dimension");
    }
  }
  //
  correlation_.resize(ncenters_);

  if(ncenters_==1){
    std::vector<double> corr(1,0.0);
    if(parseVector("CORRELATION",corr,true)){
      diagonal_ = false;
    }
    correlation_[0] = corr;
  }
  else {
    for(unsigned int i=0;i<ncenters_; i++){
      std::vector<double> corr(1,0.0);
      if(parseNumberedVector("CORRELATION",(i+1),corr)){
        diagonal_ = false;
      }
      correlation_[i] = corr;
    }
  }

  if(!diagonal_ && getDimension()!=2){
    plumed_merror(getName()+": CORRELATION is only defined for two-dimensional Gaussians for now.");
  }
  for(unsigned int i=0; i<correlation_.size(); i++){
    if(correlation_[i].size()!=1){
      plumed_merror(getName()+": only one value should be given in CORRELATION");
    }
    for(unsigned int k=0;k<correlation_[i].size(); k++){
      if(correlation_[i][k] <= -1.0 ||  correlation_[i][k] >= 1.0){
        plumed_merror(getName()+": values given in CORRELATION should be between -1.0 and 1.0" );
      }
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
  checkRead();
}


double TD_Gaussian::getValue(const std::vector<double>& argument) const {
  double value=0.0;
  if(diagonal_){
    for(unsigned int i=0;i<ncenters_;i++){
      value+=weights_[i]*GaussianDiagonal(argument, centers_[i], sigmas_[i]);
    }
  }
  else if(!diagonal_ && getDimension()==2){
    for(unsigned int i=0;i<ncenters_;i++){
      value+=weights_[i]*Gaussian2D(argument, centers_[i], sigmas_[i],correlation_[i]);
    }
  }
  return value;
}


double TD_Gaussian::GaussianDiagonal(const std::vector<double>& argument, const std::vector<double>& center, const std::vector<double>& sigma, bool normalize) const {
  double value = 1.0;
  for(unsigned int k=0; k<argument.size(); k++){
    double arg=(argument[k]-center[k])/sigma[k];
    double tmp_exp = exp(-0.5*arg*arg);
    if(normalize){tmp_exp/=(sigma[k]*sqrt(2.0*pi));}
    value*=tmp_exp;
  }
  return value;
}


double TD_Gaussian::Gaussian2D(const std::vector<double>& argument, const std::vector<double>& center, const std::vector<double>& sigma, const std::vector<double>& correlation, bool normalize) const {
  double arg1 = (argument[0]-center[0])/sigma[0];
  double arg2 = (argument[1]-center[1])/sigma[1];
  double corr = correlation[0];
  double value = (arg1*arg1 + arg2*arg2 - 2.0*corr*arg1*arg2);
  value *= -1.0 / ( 2.0*(1.0-corr*corr) );
  value = exp(value);
  if(normalize){
    value /=  2*pi*sigma[0]*sigma[1]*sqrt(1.0-corr*corr);
  }
  return value;
}

}
}
