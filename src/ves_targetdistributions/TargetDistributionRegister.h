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
#ifndef __PLUMED_ves_targetdistributions_TargetDistributionRegister_h
#define __PLUMED_ves_targetdistributions_TargetDistributionRegister_h

#include "TargetDistributionBase.h"

#include <string>
#include <cstring>
#include <vector>
#include <map>

namespace PLMD{



class TargetDistributionRegister{
private:
/// Pointer to a function which, given the type for a ReferenceConfiguration, creates it
  typedef TargetDistributionBase*(*creator_pointer)(const TargetDistributionOptions&);
/// The set of possible target distribution we can use
  std::map<std::string,creator_pointer> m;
public:
/// The destructor
  ~TargetDistributionRegister();
/// Add a new target distribution to the register
  void add( std::string type, creator_pointer );
/// Remove a target distribution to the register
  void remove(creator_pointer f);
/// Verify if a target distribution is present in the register
  bool check(std::string type);
/// Create a target distribution object
  TargetDistributionBase* create( const TargetDistributionOptions& to );
};

TargetDistributionRegister& targetDistributionRegister();

#define VARIATIONAL_REGISTER_TARGET_DISTRIBUTION(classname,type) \
  static class classname##RegisterMe{ \
    static TargetDistributionBase * create(const TargetDistributionOptions&to){return new classname(to);} \
  public: \
    classname##RegisterMe(){targetDistributionRegister().add(type,create);}; \
    ~classname##RegisterMe(){targetDistributionRegister().remove(create);}; \
  } classname##RegisterMeObject;

}
#endif
