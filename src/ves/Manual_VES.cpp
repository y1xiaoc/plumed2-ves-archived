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

#include "TargetDistributionRegister.h"

#include "cltools/CLTool.h"
#include "cltools/CLToolRegister.h"
#include "tools/Tools.h"
#include "config/Config.h"
#include <cstdio>
#include <string>
#include <vector>
#include <iostream>

namespace PLMD {
namespace ves{

//+PLUMEDOC VES_CLTOOLS manual_ves
/*
manual is a tool that you can use to construct the manual page for
a particular target distribution included in VES.

The manual constructed by this action is in html. In all probability you will never need to use this
tool. However, it is used within the scripts that generate plumed's html manual.  If you need to use this
tool outside those scripts the input is specified using the following command line arguments.

\par Examples

The following generates the html manual for the action TD_GAUSSIAN.
\verbatim
plumed manual_ves --targetdist TD_GAUSSIAN
\endverbatim


*/
//+ENDPLUMEDOC

class Manual_VES:
public CLTool
{
public:
  static void registerKeywords( Keywords& keys );
  explicit Manual_VES(const CLToolOptions& co );
  int main(FILE* in, FILE*out,Communicator& pc);
  std::string description()const{
    return "print out a description of the keywords for an action in html";
  }
};

PLUMED_REGISTER_CLTOOL(Manual_VES,"manual_ves")

void Manual_VES::registerKeywords( Keywords& keys ){
  CLTool::registerKeywords( keys );
  keys.add("compulsory","--targetdist","print the manual for this particular action");
  keys.addFlag("--vim",false,"print the keywords in vim syntax");
}

Manual_VES::Manual_VES(const CLToolOptions& co ):
CLTool(co)
{
  inputdata=commandline;
}

int Manual_VES::main(FILE* in, FILE*out,Communicator& pc){
  std::string targetdist;
  if( !parse("--targetdist",targetdist) ) return 1;
  std::cerr<<"LIST OF DOCUMENTED TARGET DISTRIBUTIONS:\n";
  std::cerr<<targetDistributionRegister()<<"\n";
  bool vimout; parseFlag("--vim",vimout);
  if( !targetDistributionRegister().printManual(targetdist,vimout)){
    fprintf(stderr,"specified action is not registered\n");
    return 1;
  }
  return 0;
}

} // End of namespace
}
