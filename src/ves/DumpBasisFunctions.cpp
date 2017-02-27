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

#include "BasisFunctions.h"
#include "TargetDistribution.h"
#include "TargetDistributionRegister.h"
#include "CoeffsVector.h"

#include "core/ActionRegister.h"
#include "core/ActionSet.h"
#include "core/PlumedMain.h"
#include "tools/File.h"
#include "tools/Grid.h"


namespace PLMD {
namespace ves {

//+PLUMEDOC VES_UTILS DUMP_BASISFUNCTIONS
/*
Dump basis functions to file.

\par Examples

*/
//+ENDPLUMEDOC


class DumpBasisFunctions :
  public Action
{
  std::vector<BasisFunctions*> bf_pntrs;
public:
  explicit DumpBasisFunctions(const ActionOptions&);
  TargetDistribution* setupTargetDistPntr(std::string keyword) const;
  void calculate() {}
  void apply() {}
  static void registerKeywords(Keywords& keys);
};


PLUMED_REGISTER_ACTION(DumpBasisFunctions,"DUMP_BASISFUNCTIONS")

void DumpBasisFunctions::registerKeywords(Keywords& keys) {
  Action::registerKeywords(keys);
  keys.add("compulsory","BASIS_SET","the label of the basis set that you want to use");
  keys.add("optional","GRID_BINS","the number of bins used for the grid for writing the basis function values and derivatives. The default value is 1000.");
  keys.add("optional","GRID_MIN","the minimum of the grid for writing the basis function values and derivatives. By default it is the minimum of the interval on which the basis functions are defined.");
  keys.add("optional","GRID_MAX","the maximum of the grid for writing the basis function values and derivatives. By default it is the maximum of the interval on which the basis functions are defined.");
  keys.add("optional","FILE_VALUES","filename of the file on which the basis function values are written. By default it is BASIS_SET.values.data.");
  keys.add("optional","FILE_DERIVS","filename of the file on which the basis function derivatives are written. By default it is BASIS_SET.derivs.data.");
  keys.add("optional","FORMAT_VALUES_DERIVS","the numerical format of the basis function values and derivatives written to file. By default it is %15.8f.\n");
  keys.add("optional","FILE_TARGETDIST_AVERAGES","filename of the file on which the averages over the target distributions are written. By default it is BASIS_SET.targetdist-averages.data.");
  keys.add("optional","FORMAT_TARGETDIST_AVERAGES","the numerical format of the target distribution averages written to file. By default it is %15.8f.\n");
  keys.add("optional","FILE_TARGETDIST","filename of the files on which the target distributions are written. By default it is BASIS_SET.targetdist-#.data.");
  keys.add("numbered","TARGET_DISTRIBUTION","the target distribution to be used.");
  keys.addFlag("IGNORE_PERIODICITY",false,"if the periodicity of the basis functions should be ignored.");
  keys.addFlag("NUMERICAL_DERIVATIES",false,"if the derivatives of the basis functions should be calculated numerically.");
}

DumpBasisFunctions::DumpBasisFunctions(const ActionOptions&ao):
  Action(ao),
  bf_pntrs(1)
{
  std::string basisset_label="";
  parse("BASIS_SET",basisset_label);
  bf_pntrs[0]=plumed.getActionSet().selectWithLabel<BasisFunctions*>(basisset_label);
  plumed_massert(bf_pntrs[0]!=NULL,"basis function "+basisset_label+" does not exist. NOTE: the basis functions should always be defined BEFORE the DUMP_BASISFUNCTIONS action.");

  unsigned int nbins = 1000;
  parse("GRID_BINS",nbins);

  std::string min_str = bf_pntrs[0]->intervalMinStr();
  std::string max_str = bf_pntrs[0]->intervalMaxStr();
  parse("GRID_MIN",min_str);
  parse("GRID_MAX",max_str);

  std::string fname_values = bf_pntrs[0]->getLabel()+".values.data";
  parse("FILE_VALUES",fname_values);
  std::string fname_derives = bf_pntrs[0]->getLabel()+".derivs.data";
  parse("FILE_DERIVS",fname_derives);
  std::string fname_targetdist_aver = bf_pntrs[0]->getLabel()+".targetdist-averages.data";
  parse("FILE_TARGETDIST_AVERAGES",fname_targetdist_aver);
  std::string fname_targetdist = bf_pntrs[0]->getLabel()+".targetdist-.data";
  parse("FILE_TARGETDIST",fname_targetdist);

  std::string fmt_values_derivs = "%15.8f";
  parse("FORMAT_VALUES_DERIVS",fmt_values_derivs);
  std::string fmt_targetdist_aver = "%15.8f";
  parse("FORMAT_TARGETDIST_AVERAGES",fmt_targetdist_aver);

  bool ignore_periodicity = false;
  parseFlag("IGNORE_PERIODICITY",ignore_periodicity);

  bool numerical_deriv = false;
  parseFlag("NUMERICAL_DERIVATIES",numerical_deriv);

  std::vector<std::string> targetdist_keywords;
  std::string str_ps="";
  for(int i=1;; i++) {
    if(!parseNumbered("TARGET_DISTRIBUTION",i,str_ps)) {break;}
    targetdist_keywords.push_back(str_ps);
  }
  checkRead();
  //
  OFile ofile_values;
  ofile_values.link(*this);
  ofile_values.enforceBackup();
  ofile_values.open(fname_values);
  OFile ofile_derivs;
  ofile_derivs.link(*this);
  ofile_derivs.enforceBackup();
  ofile_derivs.open(fname_derives);
  bf_pntrs[0]->writeBasisFunctionsToFile(ofile_values,ofile_derivs,min_str,max_str,nbins,ignore_periodicity,fmt_values_derivs,numerical_deriv);
  ofile_values.close();
  ofile_derivs.close();
  //
  std::vector<std::string> grid_min(1); grid_min[0]=bf_pntrs[0]->intervalMinStr();
  std::vector<std::string> grid_max(1); grid_max[0]=bf_pntrs[0]->intervalMaxStr();
  std::vector<unsigned int> grid_bins(1); grid_bins[0]=nbins;
  std::vector<Value*> arguments(1);
  arguments[0]= new Value(NULL,"arg",false);
  if(bf_pntrs[0]->arePeriodic() && !ignore_periodicity) {
    arguments[0]->setDomain(bf_pntrs[0]->intervalMinStr(),bf_pntrs[0]->intervalMaxStr());
  }
  else {
    arguments[0]->setNotPeriodic();
  }

  OFile ofile_targetdist_aver;
  ofile_targetdist_aver.link(*this);
  ofile_targetdist_aver.enforceBackup();
  ofile_targetdist_aver.open(fname_targetdist_aver);

  for(unsigned int i=0; i<targetdist_keywords.size(); i++) {
    std::string is; Tools::convert(i+1,is);
    //
    TargetDistribution* targetdist_pntr = setupTargetDistPntr(targetdist_keywords[i]);
    if(targetdist_pntr!=NULL) {
      targetdist_pntr->setupGrids(arguments,grid_min,grid_max,grid_bins);
      plumed_massert(targetdist_pntr->getDimension()==1,"the target distribution must be one dimensional");
      targetdist_pntr->update();
    }
    //
    std::vector<double> bf_integrals = bf_pntrs[0]->getTargetDistributionIntegrals(targetdist_pntr);
    CoeffsVector targetdist_averages = CoeffsVector("aver.targetdist-"+is,arguments,bf_pntrs,comm,false);
    targetdist_averages.setValues(bf_integrals);
    if(fmt_targetdist_aver.size()>0) {targetdist_averages.setOutputFmt(fmt_targetdist_aver);}
    targetdist_averages.writeToFile(ofile_targetdist_aver,true);
    if(targetdist_pntr!=NULL) {
      Grid* targetdist_grid_pntr = targetdist_pntr->getTargetDistGridPntr();
      std::string fname = FileBase::appendSuffix(fname_targetdist,is);
      OFile ofile;
      ofile.link(*this);
      ofile.enforceBackup();
      ofile.open(fname);
      targetdist_grid_pntr->writeToFile(ofile);
      ofile.close();
    }
    delete targetdist_pntr;
  }
  ofile_targetdist_aver.close();
  delete arguments[0]; arguments.clear();



}


TargetDistribution* DumpBasisFunctions::setupTargetDistPntr(std::string keyword) const {
  std::vector<std::string> words = Tools::getWords(keyword);
  TargetDistribution* pntr = NULL;
  if(words[0]=="DEFAULT_UNIFORM" && words.size()==1) {
    pntr = NULL;
  }
  else {
    pntr = targetDistributionRegister().create(words);
  }
  return pntr;
}







}
}
