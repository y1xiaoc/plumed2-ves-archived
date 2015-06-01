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
#ifndef __PLUMED_variational_TargetDistribution1DimRegister_h
#define __PLUMED_variational_TargetDistribution1DimRegister_h

#include <string>
#include <cstring>
#include <vector>
#include <map>
#include "TargetDistribution1DimBase.h"

namespace PLMD{


namespace variational{

class TargetDistribution1DimRegister{
private:
/// Pointer to a function which, given the type for a ReferenceConfiguration, creates it
  typedef TargetDistribution1DimBase*(*creator_pointer)(const TargetDistribution1DimOptions&);
/// The set of possible landmark selection algorithms we can work with
  std::map<std::string,creator_pointer> m;
public:
/// The destructor
  ~TargetDistribution1DimRegister();
/// Add a new landmark selection style to the register of landmark selectors
  void add( std::string type, creator_pointer );
/// Remove a landmark selection style from the register of metrics
  void remove(creator_pointer f);
/// Verify if a landmark selection style is present in the register
  bool check(std::string type);
/// Create a landmark selection object
  TargetDistribution1DimBase* create( const TargetDistribution1DimOptions& to );
};

TargetDistribution1DimRegister& targetDistribution1DimRegister();

#define VARIATIONAL_REGISTER_TARGET_DISTRIBUTION_1D(classname,type) \
  static class classname##RegisterMe{ \
    static TargetDistribution1DimBase * create(const TargetDistribution1DimOptions&to){return new classname(to);} \
  public: \
    classname##RegisterMe(){targetDistribution1DimRegister().add(type,create);}; \
    ~classname##RegisterMe(){targetDistribution1DimRegister().remove(create);}; \
  } classname##RegisterMeObject;

}
}
#endif
