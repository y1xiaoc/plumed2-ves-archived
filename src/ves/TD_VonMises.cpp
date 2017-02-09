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
#include "tools/Tools.h"
#include "GridIntegrationWeights.h"

#include <iostream>



namespace PLMD{
namespace ves{

//+PLUMEDOC VES_TARGETDIST VON_MISES
/*
Von Mises target distribution (static).

\par Examples

*/
//+ENDPLUMEDOC

class TD_VonMises: public TargetDistribution {
  // properties of the Gaussians
  std::vector< std::vector<double> > sigmas_;
  std::vector< std::vector<double> > kappas_;
  std::vector< std::vector<double> > centers_;
  std::vector< std::vector<double> > normalization_;
  std::vector<double> weights_;
  std::vector<double> periods_;
  unsigned int ncenters_;
  double VonMisesDiagonal(const std::vector<double>&, const std::vector<double>&, const std::vector<double>&, const std::vector<double>&, const std::vector<double>&) const;
  double getNormalization(const double, const double) const;
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
  keys.use("BIAS_CUTOFF");
  keys.use("WELLTEMPERED_FACTOR");
  keys.use("SHIFT_TO_ZERO");
  keys.use("NORMALIZE");
}


TD_VonMises::TD_VonMises( const TargetDistributionOptions& to ):
TargetDistribution(to),
sigmas_(0),
centers_(0),
normalization_(0),
weights_(0),
periods_(0),
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
    plumed_merror(getName()+": CENTER and SIGMA keywords seem to be missing. Note that numbered keywords start at CENTER1 and SIGMA1.");
  }
  //
  setDimension(centers_[0].size());
  ncenters_ = centers_.size();
  //
  // check centers and sigmas
  for(unsigned int i=0; i<ncenters_; i++) {
    if(centers_[i].size()!=getDimension()){plumed_merror(getName()+": one of the CENTER keyword does not match the given dimension");}
    if(sigmas_[i].size()!=getDimension()){plumed_merror(getName()+": one of the SIGMA keyword does not match the given dimension");}
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
  if(centers_.size()!=weights_.size()){plumed_merror(getName() + ": there has to be as many weights given in WEIGHTS as numbered CENTER keywords");}
  //
  if(!parseVector("PERIODS",periods_,true)){periods_.assign(getDimension(),2*pi);}
  if(periods_.size()!=getDimension()){plumed_merror(getName() + ": the number of values given in PERIODS does not match the dimension of the distribution");}
  //
  double sum_weights=0.0;
  for(unsigned int i=0;i<weights_.size();i++){sum_weights+=weights_[i];}
  for(unsigned int i=0;i<weights_.size();i++){weights_[i]/=sum_weights;}
  //
  normalization_.resize(ncenters_);
  for(unsigned int i=0; i<ncenters_; i++){
    normalization_[i].resize(getDimension());
    for(unsigned int k=0; k<getDimension(); k++){
      normalization_[i][k] = getNormalization(kappas_[i][k],periods_[k]);
    }
  }
  checkRead();
}


double TD_VonMises::getValue(const std::vector<double>& argument) const {
  double value=0.0;
  for(unsigned int i=0;i<ncenters_;i++){
    value+=weights_[i]*VonMisesDiagonal(argument, centers_[i], kappas_[i],periods_,normalization_[i]);
  }
  return value;
}


double TD_VonMises::VonMisesDiagonal(const std::vector<double>& argument, const std::vector<double>& center, const std::vector<double>& kappa, const std::vector<double>& periods, const std::vector<double>& normalization) const {
  double value = 1.0;
  for(unsigned int k=0; k<argument.size(); k++){
    double arg = kappa[k]*cos( ((2*pi)/periods[k])*(argument[k]-center[k]) );
    value*=normalization[k]*exp(arg);
  }
  return value;
}


double TD_VonMises::getNormalization(const double kappa, const double period) const {
  //
  std::vector<double> centers(1);
  centers[0] = 0.0;
  std::vector<double> kappas(1);
  kappas[0] = kappa;
  std::vector<double> periods(1);
  periods[0] = period;
  std::vector<double> norm(1);
  norm[0] = 1.0;
  //
  const unsigned int nbins = 1001;
  std::vector<double> points;
  std::vector<double> weights;
  double min = 0.0;
  double max = period;
  GridIntegrationWeights::getOneDimensionalIntegrationPointsAndWeights(points,weights,nbins,min,max);
  //
  double sum = 0.0;
  for(unsigned int l=0; l<nbins; l++){
    std::vector<double> arg(1); arg[0]= points[l];
    sum += weights[l] * VonMisesDiagonal(arg,centers,kappas,periods,norm);
  }
  return 1.0/sum;
}


}
}
