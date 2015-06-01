/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2013-2015 The plumed team
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
#include "OptimizerBase.h"
#include "OptimizerRegister.h"

namespace PLMD {

OptimizerOptions::OptimizerOptions( const std::vector<std::string>& input):
words(input)
{
}

OptimizerBase::OptimizerBase( const OptimizerOptions& to ):
type(to.words[0]),
input(to.words)
{
  input.erase( input.begin() );
}

OptimizerBase::~OptimizerBase(){
}

void OptimizerBase::parseFlag(const std::string& key, bool& t){
  Tools::parseFlag(input,key,t);
}

void OptimizerBase::checkRead() const {
  if(!input.empty()){
     std::string msg="cannot understand the following words from landmark selection input : ";
     for(unsigned i=0;i<input.size();++i) msg = msg + input[i] + ", ";
     plumed_merror(msg); 
  }
}

std::string OptimizerBase::description(){
  std::string str="Type: " + type;
  return str;
}

}

