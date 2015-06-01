/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2011-2014 The plumed team
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
#include "TargetDistribution1DimRegister.h"
#include <iostream>

namespace PLMD{

TargetDistribution1DimRegister::~TargetDistribution1DimRegister(){
  if(m.size()>0){
    std::string names="";
    for(std::map<std::string,creator_pointer>::iterator p=m.begin();p!=m.end();++p) names+=p->first+" ";
    std::cerr<<"WARNING: One dimensional target distribution "+ names +" has not been properly unregistered. This might lead to memory leak!!\n";
  }
}

TargetDistribution1DimRegister& targetDistribution1DimRegister(){
  static TargetDistribution1DimRegister ans;
  return ans;
}

void TargetDistribution1DimRegister::remove(creator_pointer f){
  for(std::map<std::string,creator_pointer>::iterator p=m.begin();p!=m.end();++p){
    if((*p).second==f){
      m.erase(p); break;
    }
  }
}

void TargetDistribution1DimRegister::add( std::string type, creator_pointer f ){
  plumed_massert(m.count(type)==0,"type has already been registered");
  m.insert(std::pair<std::string,creator_pointer>(type,f));
}

bool TargetDistribution1DimRegister::check(std::string type){
  if( m.count(type)>0 ) return true;
  return false;
}

TargetDistribution1DimBase* TargetDistribution1DimRegister::create( const TargetDistribution1DimOptions& to ){
  TargetDistribution1DimBase* lselect;
  if( check(to.words[0]) ){
     lselect=m[to.words[0]](to);
     lselect->checkRead();
  } else lselect=NULL;
  return lselect;
}

}
