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
#include "BasisFunctions.h"

#include "core/ActionRegister.h"
#include "core/ActionSet.h"
#include "core/PlumedMain.h"
#include "tools/File.h"



// using namespace std;

namespace PLMD{
namespace generic{

//+PLUMEDOC FUNCTION DUMP_BASISFUNCTIONS
/*

*/
//+ENDPLUMEDOC


class DumpBasisFunctions :
  public Action
{
  BasisFunctions* bf_pointer;
public:
  explicit DumpBasisFunctions(const ActionOptions&);
  void calculate(){}
  void apply(){}
  static void registerKeywords(Keywords& keys);
};


PLUMED_REGISTER_ACTION(DumpBasisFunctions,"DUMP_BASISFUNCTIONS")

void DumpBasisFunctions::registerKeywords(Keywords& keys){
  Action::registerKeywords(keys);
  keys.add("compulsory","BASIS_SET","the label of the basis set that you want to use");
  keys.add("optional","GRID_BINS","the number of bins used for the grid. The default value is 1000.");
  keys.add("optional","FILE_VALUES","filename of the file on which the basis function values are written. By default it is LABEL.values.data.");
  keys.add("optional","FILE_DERIVS","filename of the file on which the basis function derivatives are written. By default it is LABEL.derivs.data.");
}

DumpBasisFunctions::DumpBasisFunctions(const ActionOptions&ao):
Action(ao)
{
  std::string basisset_label="";
  parse("BASIS_SET",basisset_label);
  bf_pointer=plumed.getActionSet().selectWithLabel<BasisFunctions*>(basisset_label);

  unsigned int nbins = 1000;
  parse("GRID_BINS",nbins);

  std::string fname_values = bf_pointer->getLabel()+".values.data";
  parse("FILE_VALUES",fname_values);
  std::string fname_derives = bf_pointer->getLabel()+".derivs.data";
  parse("FILE_DERIVS",fname_derives);
  checkRead();

  OFile ofile_values;
  ofile_values.link(*this);
  ofile_values.open(fname_values);
  OFile ofile_derivs;
  ofile_derivs.link(*this);
  ofile_derivs.open(fname_derives);

  bf_pointer->writeBasisFunctionsToFile(ofile_values,ofile_derivs,nbins);

  ofile_values.close();
  ofile_derivs.close();
}


}
}
