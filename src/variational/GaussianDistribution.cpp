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

class GaussianDistribution: public TargetDistributionBase {
  // properties of the Gaussians
  std::vector< std::vector<double> > sigmas;
  std::vector< std::vector<double> > centers;
  std::vector<double> weights;
  bool normalize_distribution;
  unsigned int ngaussians;
  double GaussianDist(const std::vector<double>, const std::vector<double>, const std::vector<double>, const bool normalize=true);
public:
  static void registerKeywords(Keywords&);
  GaussianDistribution(const TargetDistributionOptions& to);
  double distribution(const std::vector<double>);
};


VARIATIONAL_REGISTER_TARGET_DISTRIBUTION(GaussianDistribution,"GAUSSIAN")


void GaussianDistribution::registerKeywords(Keywords& keys){
  TargetDistributionBase::registerKeywords(keys);
  keys.add("numbered","CENTER","The centers of the Gaussians.");
  keys.add("numbered","SIGMA","The sigmas of the Gaussians.");
  keys.add("optional","WEIGHTS","The weights of the Gaussians.");
  keys.addFlag("DO_NOT_NORMALIZE",false,"If the distribution should not be normalized.");
  keys.remove("DIMENSION");
}


GaussianDistribution::GaussianDistribution( const TargetDistributionOptions& to ):
TargetDistributionBase(to)
{
  for(unsigned int i=0;; i++) {
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
  // check centers and sigmas
  for(unsigned int i=0; i<ngaussians; i++) {
    plumed_massert(centers[i].size()==getDimension(),"one of the CENTER keyword does not match the given dimension");
    plumed_massert(sigmas[i].size()==getDimension(),"one of the CENTER keyword does not match the given dimension");
  }
  //
  if(!parseVector("WEIGHTS",weights,true)){weights.assign(centers.size(),1.0);}
  plumed_massert(centers.size()==weights.size(),"there has to be as many weights given in WEIGHTS as numbered CENTER keywords");
  //
  bool do_not_normalize;
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


double GaussianDistribution::distribution(const std::vector<double> argument) {
  double value=0.0;
  for(unsigned int i=0;i<ngaussians;i++){
    value+=weights[i]*GaussianDist(argument, centers[i], sigmas[i],normalize_distribution);
  }
  return value;
}


double GaussianDistribution::GaussianDist(const std::vector<double> argument, const std::vector<double> center, const std::vector<double> sigma, bool normalize){
  double value = 1.0;
  for(unsigned int k=0; k<argument.size(); k++){
    double argT=(argument[k]-center[k])/sigma[k];
    double tmp_exp = exp(-0.5*argT*argT);
    if(normalize){tmp_exp/=(sigma[k]*sqrt(2.0*pi));}
    value*=tmp_exp;
  }
  return value;
}


}
