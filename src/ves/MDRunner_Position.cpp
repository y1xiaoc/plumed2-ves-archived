/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2013 The plumed team
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

#include "colvar/Colvar.h"
#include "core/ActionRegister.h"

#include <string>
#include <cmath>


namespace PLMD{
namespace ves{

class MDRunner_Position : public Colvar {
public:
  static void registerKeywords( Keywords& keys );
  explicit MDRunner_Position( const ActionOptions& );
  virtual void calculate();
};

//+PLUMEDOC VES_COLVAR MDRUNNER_POSITION
/*
Use as CVs the postions from dynamics of \ref md_linearexpansion

\par Examples

*/
//+ENDPLUMEDOC

PLUMED_REGISTER_ACTION(MDRunner_Position,"MDRUNNER_POSITION")

void MDRunner_Position::registerKeywords( Keywords& keys ){
  Colvar::registerKeywords( keys );
  keys.addOutputComponent("x","default","x-postion of the walker");
  keys.addOutputComponent("y","default","y-postion of the walker");
  keys.addOutputComponent("z","default","z-postion of the walker");
}

MDRunner_Position::MDRunner_Position(const ActionOptions& ao):
PLUMED_COLVAR_INIT(ao)
{
  checkRead();
  log.printf("  using x,y,z coordinates of a MD Runner\n");

  addComponentWithDerivatives("x"); componentIsNotPeriodic("x");
  addComponentWithDerivatives("y"); componentIsNotPeriodic("y");
  addComponentWithDerivatives("z"); componentIsNotPeriodic("z");
  std::vector<AtomNumber> atoms(1); atoms[0].setSerial(1);
  requestAtoms(atoms);
}

void MDRunner_Position::calculate(){
  Vector pos = getPosition(0);
  Value* valuex=getPntrToComponent("x");
  Value* valuey=getPntrToComponent("y");
  Value* valuez=getPntrToComponent("z");

  valuex->set(pos[0]); setAtomsDerivatives(valuex,0,Vector(1,0,0));
  valuey->set(pos[1]); setAtomsDerivatives(valuey,0,Vector(0,1,0));
  valuez->set(pos[2]); setAtomsDerivatives(valuez,0,Vector(0,0,1));
}

}
}
