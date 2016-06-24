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
#include "TargetDistributionRegister.h"

#include <iostream>

namespace PLMD{

TargetDistributionRegister::~TargetDistributionRegister(){
  if(m.size()>0){
    std::string names="";
    for(std::map<std::string,creator_pointer>::iterator p=m.begin();p!=m.end();++p) names+=p->first+" ";
    std::cerr<<"WARNING: One dimensional target distribution "+ names +" has not been properly unregistered. This might lead to memory leak!!\n";
  }
}

TargetDistributionRegister& targetDistributionRegister(){
  static TargetDistributionRegister ans;
  return ans;
}

void TargetDistributionRegister::remove(creator_pointer f){
  for(std::map<std::string,creator_pointer>::iterator p=m.begin();p!=m.end();++p){
    if((*p).second==f){
      m.erase(p); break;
    }
  }
}

void TargetDistributionRegister::add( std::string type, creator_pointer f ){
  plumed_massert(m.count(type)==0,"type has already been registered");
  m.insert(std::pair<std::string,creator_pointer>(type,f));
}

bool TargetDistributionRegister::check(std::string type){
  if( m.count(type)>0 ) return true;
  return false;
}

TargetDistribution* TargetDistributionRegister::create( const TargetDistributionOptions& to ){
  TargetDistribution* lselect;
  if( check(to.words[0]) ){
     lselect=m[to.words[0]](to);
     lselect->checkRead();
  }
  else{
    plumed_merror("problem with setting up a target distribution, distribution of the type " + to.words[0] + " does not exist");
  }
  return lselect;
}

}
