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

#include "CounterBase.h"
#include "tools/File.h"

namespace PLMD{

CounterBase::CounterBase(const bool counter_active):
counter(0),
field_name_("iteration"),
isActive(counter_active) {
}


void CounterBase::turnOnCounter() {
  isActive=true;
}


void CounterBase::turnOffCounter() {
  isActive=false;
}


bool  CounterBase::isCounterActive() {
  return isActive;
}


void CounterBase::resetCounter() {
  counter=0;
}


void CounterBase::increaseCounter() {
  counter=+1;
}


void CounterBase::addToCounter(const unsigned int value) {
  counter=+value;
}


void CounterBase::setCounter(const unsigned int value) {
  counter=value;
}


unsigned int CounterBase::getCounter() const {
  return counter;
}


void CounterBase::setFieldName(const std::string field_name) {
  field_name_=field_name;
}


std::string CounterBase::getFieldName() const {
  return field_name_;
}


bool CounterBase::getCounterFieldFromFile(IFile& ifile) {
  bool field_found=false;
  if(ifile.FieldExist(field_name_)){
    field_found=true;
    int int_tmp;
    ifile.scanField(field_name_,int_tmp);
    counter=(unsigned int) int_tmp;
  }
  return field_found;
}


bool CounterBase::isCounterFieldInFile(IFile& ifile) {
  return ifile.FieldExist(field_name_);
}


void CounterBase::writeCounterFieldToFile(OFile& ofile) {
  if(isActive){
    ofile.addConstantField(field_name_).printField(field_name_,(int) counter);
  }
}


}
