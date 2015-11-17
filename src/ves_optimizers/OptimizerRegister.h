/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2013,2014 The plumed team
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
#ifndef __PLUMED_variational_OptimizerRegister_h
#define __PLUMED_variational_OptimizerRegister_h

#include "OptimizerBase.h"

#include <string>
#include <cstring>
#include <vector>
#include <map>

namespace PLMD{



class OptimizerRegister{
private:
/// Pointer to a function which, given the type for a ReferenceConfiguration, creates it
  typedef OptimizerBase*(*creator_pointer)(const OptimizerOptions&);
/// The set of possible optimization algorithms we can use
  std::map<std::string,creator_pointer> m;
public:
/// The destructor
  ~OptimizerRegister();
/// Add a new optimization algorithm to the register
  void add( std::string type, creator_pointer );
/// Remove a optimization algorithm to the register
  void remove(creator_pointer f);
/// Verify if a algorithm is present
  bool check(std::string type);
/// Create a optimizer object
  OptimizerBase* create( const OptimizerOptions& to );
};

OptimizerRegister& optimizerRegister();

#define VARIATIONAL_REGISTER_OPTIMIZER(classname,type) \
  static class classname##RegisterMe{ \
    static OptimizerBase * create(const OptimizerOptions&to){return new classname(to);} \
  public: \
    classname##RegisterMe(){optimizerRegisterRegister().add(type,create);}; \
    ~classname##RegisterMe(){optimizerRegisterRegister().remove(create);}; \
  } classname##RegisterMeObject;

}
#endif
