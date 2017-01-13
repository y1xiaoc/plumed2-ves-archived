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
#ifndef __PLUMED_ves_TargetDistributionRegister_h
#define __PLUMED_ves_TargetDistributionRegister_h

#include "TargetDistribution.h"

#include <string>
#include <cstring>
#include <vector>
#include <map>


namespace PLMD{
namespace ves{



class TargetDistributionRegister{
private:
/// Pointer to a function which, given the type for a TargetDistribution, creates it
  typedef TargetDistribution*(*creator_pointer)(const TargetDistributionOptions&);
  //
  typedef void(*keywords_pointer)(Keywords&);
/// The set of possible target distribution we can use
  std::map<std::string,creator_pointer> m;
  //
  std::map<std::string,Keywords> mk;
  //
  typedef std::map<std::string,creator_pointer>::const_iterator const_mIterator;
public:
/// The destructor
  ~TargetDistributionRegister();
/// Add a new target distribution to the register
  void add( std::string type, creator_pointer, keywords_pointer kf);
/// Remove a target distribution to the register
  void remove(creator_pointer f);
/// Verify if a target distribution is present in the register
  bool check(std::string type);
/// Create a target distribution object
  TargetDistribution* create( const TargetDistributionOptions& to );
  //
  std::vector<std::string> list() const;
  //
  bool printManual(const std::string& type, const bool& vimout);
};

TargetDistributionRegister& targetDistributionRegister();

std::ostream & operator<<(std::ostream &log,const TargetDistributionRegister&ar);

#define VES_REGISTER_TARGET_DISTRIBUTION(classname,type) \
  static class classname##RegisterMe{ \
    static TargetDistribution * create(const TargetDistributionOptions&to){return new classname(to);} \
  public: \
    classname##RegisterMe(){targetDistributionRegister().add(type,create,classname::registerKeywords);}; \
    ~classname##RegisterMe(){targetDistributionRegister().remove(create);}; \
  } classname##RegisterMeObject;

}
}

#endif
