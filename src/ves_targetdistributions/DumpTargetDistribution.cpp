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

#include "core/ActionRegister.h"
#include "core/ActionSet.h"
#include "core/PlumedMain.h"
#include "tools/File.h"
#include "tools/Grid.h"




// using namespace std;

namespace PLMD{
namespace generic{

//+PLUMEDOC FUNCTION DUMP_BASISFUNCTIONS
/*

*/
//+ENDPLUMEDOC


class DumpTargetDistribution :
  public Action
{
public:
  explicit DumpTargetDistribution(const ActionOptions&);
  void calculate(){}
  void apply(){}
  static void registerKeywords(Keywords& keys);
};


PLUMED_REGISTER_ACTION(DumpTargetDistribution,"DUMP_TARGET_DISTRIBUTION")

void DumpTargetDistribution::registerKeywords(Keywords& keys){
  Action::registerKeywords(keys);
  keys.add("compulsory","GRID_MIN","the lower bounds for the grid");
  keys.add("compulsory","GRID_MAX","the upper bounds for the grid");
  keys.add("compulsory","GRID_BINS","the number of bins used for the grid.");
  keys.add("optional","GRID_PERIODICITY","specfiy if the individual arguments should be made periodic (YES) or not (NO). By default all arguments are not periodic.");
  keys.add("compulsory","FILE","filename of the files on which the target distribution are written.");
  keys.add("compulsory","TARGET_DISTRIBUTION","the target distribution to be used.");
}

DumpTargetDistribution::DumpTargetDistribution(const ActionOptions&ao):
Action(ao)
{
  std::string targetdist_keyword;
  parse("TARGET_DISTRIBUTION",targetdist_keyword);
  std::vector<std::string> words = Tools::getWords(targetdist_keyword);
  TargetDistribution* targetdist_pntr = targetDistributionRegister().create(words);

  std::string fname;
  parse("FILE",fname);

  unsigned int nargs = targetdist_pntr->getDimension();

  std::vector<unsigned int> grid_bins(nargs);
  parseVector("GRID_BINS",grid_bins);
  plumed_massert(grid_bins.size()==nargs,"mismatch between number of values given in GRID_BINS and dimension of target distribution");
  std::vector<std::string> grid_min(nargs);
  parseVector("GRID_MIN",grid_min);
  plumed_massert(grid_min.size()==nargs,"mismatch between number of values given in GRID_MIN and dimension of target distribution");
  std::vector<std::string> grid_max(nargs);
  parseVector("GRID_MAX",grid_max);
  plumed_massert(grid_max.size()==nargs,"mismatch between number of values given in GRID_MIN and dimension of target distribution");
  std::vector<std::string> grid_periodicity(nargs);
  parseVector("GRID_PERIODICITY",grid_periodicity);
  if(grid_periodicity.size()==0){grid_periodicity.assign(nargs,"NO");}
  plumed_massert(grid_periodicity.size()==nargs,"mismatch between number of values given in GRID_PERIODICITY and dimension of target distribution");
  checkRead();
  //
  std::vector<Value*> arguments(nargs);
  for(unsigned int i=0; i < nargs; i++) {
    std::string is; Tools::convert(i+1,is);
    arguments[i]= new Value(NULL,"arg"+is,false);
    if(grid_periodicity[i]=="YES"){
      arguments[i]->setDomain(grid_min[i],grid_max[i]);
    }
    else if(grid_periodicity[i]=="NO"){
      arguments[i]->setNotPeriodic();
    }
    else{
      plumed_merror("wrong value given in GRID_PERIODICITY, either give YES or NO");
    }
  }
  //
  Grid ps_grid = Grid("targetdist",arguments,grid_min,grid_max,grid_bins,false,false);
  targetdist_pntr->calculateDistributionOnGrid(&ps_grid);

  std::vector<double> integration_weights = GridIntegrationWeights::getTrapezoidalIntegrationWeights(&ps_grid,"weights_grid.data");
  double sum_grid=0.0;
  double sum_grid2=0.0;
  for(unsigned int i=0; i<ps_grid.getSize(); i++){
    sum_grid  += ps_grid.getValue(i);
    sum_grid2 += integration_weights[i]*ps_grid.getValue(i);
  }
  sum_grid *= ps_grid.getBinVolume();
  log.printf("  target distribtion summed over the grid: %16.12f (normal sum)\n",sum_grid);
  log.printf("  target distribtion summed over the grid: %16.12f (with weights)\n",sum_grid2);
  //
  OFile ofile;
  ofile.link(*this);
  ofile.enforceBackup();
  ofile.open(fname);
  ps_grid.writeToFile(ofile);
  ofile.close();
  //
  delete targetdist_pntr;
  for(unsigned int i=0; i < nargs; i++) {
    delete arguments[i];
  }
  arguments.clear();


}





}
}
