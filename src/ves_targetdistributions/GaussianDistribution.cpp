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
#include "TargetDistribution.h"
#include "TargetDistributionRegister.h"

#include "tools/Keywords.h"

namespace PLMD {

//+PLUMEDOC INTERNAL GAUSSIAN
/*
  Gaussian target distribution
*/
//+ENDPLUMEDOC

class GaussianDistribution: public TargetDistribution {
  // properties of the Gaussians
  std::vector< std::vector<double> > sigmas;
  std::vector< std::vector<double> > centers;
  std::vector< std::vector<double> > correlation;
  std::vector<double> weights;
  bool normalize_distribution;
  bool diagonal;
  unsigned int ngaussians;
  double GaussianDiagonal(const std::vector<double>&, const std::vector<double>&, const std::vector<double>&, const bool normalize=true) const;
  double Gaussian2D(const std::vector<double>&, const std::vector<double>&, const std::vector<double>&, const std::vector<double>&, const bool normalize=true) const;
public:
  static void registerKeywords(Keywords&);
  explicit GaussianDistribution(const TargetDistributionOptions& to);
  double getValue(const std::vector<double>&) const;
};


VARIATIONAL_REGISTER_TARGET_DISTRIBUTION(GaussianDistribution,"GAUSSIAN")


void GaussianDistribution::registerKeywords(Keywords& keys){
  TargetDistribution::registerKeywords(keys);
  keys.add("numbered","CENTER","The centers of the Gaussians.");
  keys.add("numbered","SIGMA","The sigmas of the Gaussians.");
  keys.add("numbered","CORRELATION","The correlation between the arguments, currently only works for two-dimensional Gaussians ");
  keys.add("optional","WEIGHTS","The weights of the Gaussians.");
  keys.addFlag("DO_NOT_NORMALIZE",false,"If the distribution should not be normalized.");
}


GaussianDistribution::GaussianDistribution( const TargetDistributionOptions& to ):
TargetDistribution(to),
diagonal(true)
{
  for(unsigned int i=1;; i++) {
    std::vector<double> tmp_center;
    if(!parseNumberedVector("CENTER",i,tmp_center) ){break;}
    centers.push_back(tmp_center);
  }
  for(unsigned int i=0;; i++) {
    std::vector<double> tmp_sigma;
    if(!parseNumberedVector("SIGMA",i,tmp_sigma) ){break;}
    sigmas.push_back(tmp_sigma);
  }
  plumed_massert(centers.size()==sigmas.size(),"there has to be an equal amount of numbered CENTER and SIGMA keywords");
  setDimension(centers[0].size());
  ngaussians = centers.size();

  correlation.resize(ngaussians);
  //
  for(unsigned int i=0;i<ngaussians; i++){
    std::vector<double> corr;
    if(parseNumberedVector("CORRELATION",i,corr,true)){
      plumed_massert(getDimension()==2,"CORRELATION is only defined for two-dimensional Gaussians");
      plumed_massert(corr.size()==1,"only one value should be given in CORRELATION");
      for(unsigned int k=0;k<corr.size(); k++){
        plumed_massert(corr[k] >= -1.0 && corr[k] <= 1.0,"values given in CORRELATION should be between -1.0 and 1.0" );
      }
      correlation[i] = corr;
      diagonal = false;
    }
    else {
      corr.assign(1,0.0);
      correlation[i] = corr;
    }
  }
  // check centers and sigmas
  for(unsigned int i=0; i<ngaussians; i++) {
    plumed_massert(centers[i].size()==getDimension(),"one of the CENTER keyword does not match the given dimension");
    plumed_massert(sigmas[i].size()==getDimension(),"one of the CENTER keyword does not match the given dimension");
  }
  //
  if(!parseVector("WEIGHTS",weights,true)){weights.assign(centers.size(),1.0);}
  plumed_massert(centers.size()==weights.size(),"there has to be as many weights given in WEIGHTS as numbered CENTER keywords");
  //
  bool do_not_normalize=false;
  parseFlag("DO_NOT_NORMALIZE",do_not_normalize);
  normalize_distribution=!do_not_normalize;
  if(normalize_distribution){
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


double GaussianDistribution::getValue(const std::vector<double>& argument) const {
  double value=0.0;
  if(diagonal){
    for(unsigned int i=0;i<ngaussians;i++){
      value+=weights[i]*GaussianDiagonal(argument, centers[i], sigmas[i],normalize_distribution);
    }
  }
  else if(!diagonal && getDimension()==2){
    for(unsigned int i=0;i<ngaussians;i++){
      value+=weights[i]*Gaussian2D(argument, centers[i], sigmas[i],correlation[i],normalize_distribution);
    }
  }
  return value;
}


double GaussianDistribution::GaussianDiagonal(const std::vector<double>& argument, const std::vector<double>& center, const std::vector<double>& sigma, bool normalize) const {
  double value = 1.0;
  for(unsigned int k=0; k<argument.size(); k++){
    double arg=(argument[k]-center[k])/sigma[k];
    double tmp_exp = exp(-0.5*arg*arg);
    if(normalize){tmp_exp/=(sigma[k]*sqrt(2.0*pi));}
    value*=tmp_exp;
  }
  return value;
}


double GaussianDistribution::Gaussian2D(const std::vector<double>& argument, const std::vector<double>& center, const std::vector<double>& sigma, const std::vector<double>& correlation, bool normalize) const {
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
