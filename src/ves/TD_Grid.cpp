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
#include "GridIntegrationWeights.h"
#include "VesTools.h"

#include "tools/Keywords.h"
#include "tools/Grid.h"
#include "core/Value.h"
#include "tools/File.h"


namespace PLMD{
namespace ves{

class TD_Grid : public TargetDistribution {
  Grid* distGrid_;
  std::vector<double> minima_;
  std::vector<double> maxima_;
  std::vector<bool> periodic_;
  bool zero_outside_;
public:
  static void registerKeywords( Keywords&);
  explicit TD_Grid( const TargetDistributionOptions& to );
  double getValue(const std::vector<double>&) const ;
};


VES_REGISTER_TARGET_DISTRIBUTION(TD_Grid,"GRID_DIST")


void TD_Grid::registerKeywords(Keywords& keys) {
  TargetDistribution::registerKeywords(keys);
  keys.add("compulsory","FILE","the name of the grid file contaning the target distribution");
  // keys.addFlag("NOSPLINE",false,"specifies that no spline interpolation is to be used when calculating the target distribution");
  keys.addFlag("ZERO_OUTSIDE",false,"by default the target distribution is continuous such that values outside the given grid are the same as at the boundary. This can be changed by using this flag which will make values outside the grid to be taken as zero.");
}


TD_Grid::TD_Grid(const TargetDistributionOptions& to):
TargetDistribution(to),
distGrid_(NULL),
minima_(0),
maxima_(0),
zero_outside_(false)
{
  std::string filename;
  parse("FILE",filename);
  parseFlag("ZERO_OUTSIDE",zero_outside_);
  bool no_spline=false;
  // parseFlag("NOSPLINE",no_spline);
  bool use_spline = !no_spline;

  checkRead();

  std::string gridlabel;
  std::vector<std::string> arglabels;
  std::vector<std::string> argmin;
  std::vector<std::string> argmax;
  std::vector<bool> argperiodic;
  std::vector<unsigned int> argnbins;
  bool has_deriv = false;
  unsigned int nargs = VesTools::getGridFileInfo(filename,gridlabel,arglabels,argmin,argmax,argperiodic,argnbins,has_deriv);
  if(nargs==0){
    plumed_merror("Target distribution of type GRID: problem in parsing information from grid file");
  }

  setDimension(arglabels.size());
  std::vector<Value*> arguments(arglabels.size());
  for(unsigned int i=0; i < arglabels.size(); i++) {
    arguments[i]= new Value(NULL,arglabels[i],false);
    if(argperiodic[i]){
      arguments[i]->setDomain(argmin[i],argmax[i]);
    }
    else {
      arguments[i]->setNotPeriodic();
    }
  }

  IFile gridfile; gridfile.open(filename);
  if(has_deriv){
    distGrid_=Grid::create(gridlabel,arguments,gridfile,false,use_spline,true);
  }
  else {
    distGrid_=Grid::create(gridlabel,arguments,gridfile,false,false,false);
    if(use_spline){distGrid_->enableSpline();}
  }

  plumed_massert(distGrid_->getDimension()==getDimension(),"Target distribution of type GRID: mismatch in the dimension of the read-in grid and tha arguments given in ARGS");

  minima_.resize(getDimension());
  maxima_.resize(getDimension());
  periodic_.resize(getDimension());
  for (unsigned int i=0; i < getDimension(); i++) {
    Tools::convert(distGrid_->getMin()[i],minima_[i]);
    Tools::convert(distGrid_->getMax()[i],maxima_[i]);
    periodic_[i] = argperiodic[i];
    if(periodic_[i]){maxima_[i]-=distGrid_->getDx()[i];}
  }

}


double TD_Grid::getValue(const std::vector<double>& argument) const {
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
}
