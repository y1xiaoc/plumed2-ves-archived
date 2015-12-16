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
#include "tools/Grid.h"
#include "core/Value.h"
#include "tools/File.h"

namespace PLMD {

class DistributionFromGrid : public TargetDistributionBase {
  double normalization;
  Grid* distGrid;
  std::vector<double> minima;
  std::vector<double> maxima;
public:
  static void registerKeywords( Keywords&);
  explicit DistributionFromGrid( const TargetDistributionOptions& to );
  double getValue(const std::vector<double>&) const ;
};


VARIATIONAL_REGISTER_TARGET_DISTRIBUTION(DistributionFromGrid,"GRID")


void DistributionFromGrid::registerKeywords(Keywords& keys) {
  TargetDistributionBase::registerKeywords(keys);
  keys.add("compulsory","FILE","the name of the file contaning the target distribtion");
  keys.add("compulsory","ARGS","the arguments given in the grid file");
  keys.add("compulsory","LABEL","the label given in the grid file");
  // keys.addFlag("NOSPLINE",false,"specifies that no spline interpolation is to be used when calculating the target distribution");
  keys.addFlag("NORMALIZE",false,"specifies that the target distribtion should be normalized by integrating over it. Otherwise it is assumed that it is normalized.");
}


DistributionFromGrid::DistributionFromGrid(const TargetDistributionOptions& to):
TargetDistributionBase(to)
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
  distGrid=Grid::create(gridlabel,arguments,gridfile,sparsegrid,spline,false);
  plumed_massert(distGrid->getDimension()==getDimension(),"mismatch in the dimension of the read-in grid and tha arguments given in ARGS");

  minima.resize(getDimension());
  maxima.resize(getDimension());
  for (unsigned int i=0; i < getDimension(); i++) {
    Tools::convert(distGrid->getMin()[i],minima[i]);
    Tools::convert(distGrid->getMax()[i],maxima[i]);
  }
  if(normalize){
    normalization = 0.0;
    for(unsigned int l=0; l<distGrid->getSize(); l++)
    {
     normalization += distGrid->getValue(l);
    }
    normalization *= distGrid->getBinVolume();
    distGrid->scaleAllValuesAndDerivatives(1.0/normalization);
  }
  setNormalized();
}


double DistributionFromGrid::getValue(const std::vector<double>& argument) const {
  double outside = 0.0;
  for(unsigned int k=0; k<getDimension(); k++){
    if(argument[k] < minima[k] || argument[k] > maxima[k]){return outside;}
  }
  return distGrid->getValue(argument);
}


}
