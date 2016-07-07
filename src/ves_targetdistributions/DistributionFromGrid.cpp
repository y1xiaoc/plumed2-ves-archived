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
#include "ves_tools/GridIntegrationWeights.h"


#include "tools/Keywords.h"
#include "tools/Grid.h"
#include "core/Value.h"
#include "tools/File.h"

namespace PLMD {

class DistributionFromGrid : public TargetDistribution {
  double normalization_;
  Grid* distGrid_;
  std::vector<double> minima_;
  std::vector<double> maxima_;
  bool zero_outside_;
public:
  static void registerKeywords( Keywords&);
  explicit DistributionFromGrid( const TargetDistributionOptions& to );
  double getValue(const std::vector<double>&) const ;
};


VES_REGISTER_TARGET_DISTRIBUTION(DistributionFromGrid,"GRID")


void DistributionFromGrid::registerKeywords(Keywords& keys) {
  TargetDistribution::registerKeywords(keys);
  keys.add("compulsory","FILE","the name of the grid file contaning the target distribution");
  keys.add("compulsory","ARGS","the arguments given in the grid file");
  keys.add("compulsory","LABEL","the label given in the grid file");
  // keys.addFlag("NOSPLINE",false,"specifies that no spline interpolation is to be used when calculating the target distribution");
  keys.addFlag("NORMALIZE",false,"specifies that the target distribution should be normalized by integrating over it. Otherwise it is assumed that it is normalized.");
  keys.addFlag("ZERO_OUTSIDE",false,"by default the target distribution is continuous such that values outside the given grid are the same as at the boundary. This can be changed by using this flag which will make values outside the grid to be taken as zero.");
}


DistributionFromGrid::DistributionFromGrid(const TargetDistributionOptions& to):
TargetDistribution(to),
normalization_(0.0),
distGrid_(NULL),
minima_(0),
maxima_(0),
zero_outside_(false)
{
  std::string filename;
  parse("FILE",filename);
  std::string gridlabel;
  parse("LABEL",gridlabel);
  std::vector<std::string> arglabels;
  parseVector("ARGS",arglabels);
  setDimension(arglabels.size());
  bool normalize=false;
  parseFlag("NORMALIZE",normalize);
  parseFlag("ZERO_OUTSIDE",zero_outside_);
  // bool nospline=false;
  // parseFlag("NOSPLINE",nospline);
  // bool spline=!nospline;
  bool spline = false;
  bool sparsegrid=false;
  checkRead();

  std::vector<Value*> arguments(arglabels.size());
  for (unsigned int i=0; i < arglabels.size(); i++) {
    arguments[i]= new Value(NULL,arglabels[i],false);
    arguments[i]->setNotPeriodic();
  }
  IFile gridfile; gridfile.open(filename);
  distGrid_=Grid::create(gridlabel,arguments,gridfile,sparsegrid,spline,false);
  plumed_massert(distGrid_->getDimension()==getDimension(),"mismatch in the dimension of the read-in grid and tha arguments given in ARGS");

  minima_.resize(getDimension());
  maxima_.resize(getDimension());
  for (unsigned int i=0; i < getDimension(); i++) {
    Tools::convert(distGrid_->getMin()[i],minima_[i]);
    Tools::convert(distGrid_->getMax()[i],maxima_[i]);
  }

  normalization_ = 0.0;
  std::vector<double> integration_weights = GridIntegrationWeights::getIntegrationWeights(distGrid_);
  for(unsigned int l=0; l<distGrid_->getSize(); l++){
   normalization_ += integration_weights[l]*distGrid_->getValue(l);
  }

  if(normalize){
    distGrid_->scaleAllValuesAndDerivatives(1.0/normalization_);
  }

   setNormalized();
}


double DistributionFromGrid::getValue(const std::vector<double>& argument) const {
  double outside = 0.0;
  std::vector<double> arg = argument;
  for(unsigned int k=0; k<getDimension(); k++){
    if(zero_outside_ && (argument[k] < minima_[k] || argument[k] > maxima_[k])){
      return outside;
    }
    else if(argument[k] < minima_[k]){
      arg[k] = minima_[k];
    }
    else if(argument[k] > maxima_[k]){
      arg[k] =maxima_[k];
    }
  }
  return distGrid_->getValue(arg);
}


}
