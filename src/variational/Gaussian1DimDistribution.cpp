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
#include "TargetDistribution1DimBase.h"
#include "TargetDistribution1DimRegister.h"

namespace PLMD {

class Gaussian1DimDistribution : public TargetDistribution1DimBase {
 // properties of the Gaussians 
 std::vector<double> sigmas;
 std::vector<double> centers;
 std::vector<double> weights;
 bool normalize_distribution;
 double GaussianDist(const double, const double, const double, const bool normalize=true);
public:
  Gaussian1DimDistribution( const TargetDistribution1DimOptions& to );
  double distribution(const double);
};

VARIATIONAL_REGISTER_TARGET_DISTRIBUTION_1D(Gaussian1DimDistribution,"GAUSSIAN")

Gaussian1DimDistribution::Gaussian1DimDistribution( const TargetDistribution1DimOptions& to ):
TargetDistribution1DimBase(to)
{
 parseVector("CENTER",centers);
 parseVector("SIGMA",sigmas);
 plumed_massert(centers.size()==sigmas.size(),"CENTER and SIGMA need to be of the same size");
 // 
 if(!parseVector("WEIGHT",weights,true)){weights.assign(centers.size(),1.0);}
 plumed_massert(centers.size()==weights.size(),"CENTER and WEIGHT need to be of the same size");
 //
 bool do_not_normalize=false;
 parseFlag("DO_NOT_NORMALIZE",do_not_normalize);
 normalize_distribution=!do_not_normalize;
 if(normalize_distribution){
  double sum_weights=0.0;
  for(unsigned int i=0;i<weights.size();i++){sum_weights+=weights[i];}
  for(unsigned int i=0;i<weights.size();i++){weights[i]/=sum_weights;}
  setNormalized();
 }else{
  setNotNormalized();
 }
}

double Gaussian1DimDistribution::distribution(const double argument)
{
 double value=0.0;
 for(unsigned int i=0;i<weights.size();i++){
  value+=weights[i]*GaussianDist(argument, centers[i], sigmas[i],normalize_distribution);
 }
 return value;
}

double Gaussian1DimDistribution::GaussianDist(const double argument, const double center, const double sigma, bool normalize)
{
 double argumentT=(argument-center)/sigma;
 double value=exp(-0.5*argumentT*argumentT);
 if(normalize){value/=(sigma*sqrt(2.0*pi));}
 return value;
}




}
