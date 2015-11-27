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

CounterBase::CounterBase(const bool active):
counter(0),
field_name_counter_("iteration"),
field_name_time_("time_"),
isActive(active) {
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
  counter++;
}


void CounterBase::addToCounter(const unsigned int value) {
  counter+=value;
}


void CounterBase::setCounter(const unsigned int value) {
  counter=value;
}


unsigned int CounterBase::getCounter() const {
  return counter;
}


double CounterBase::getCounterDbl() const {
  return (double) counter;
}


void CounterBase::setCounterFieldName(const std::string field_name_counter) {
  field_name_counter_=field_name_counter;
}


std::string CounterBase::getCounterFieldName() const {
  return field_name_counter_;
}


void CounterBase::setTimeFieldName(const std::string field_name_time) {
  field_name_time_=field_name_time;
}


std::string CounterBase::getTimeFieldName() const {
  return field_name_time_;
}



bool CounterBase::getCounterFieldFromFile(IFile& ifile) {
  bool field_found=false;
  if(ifile.FieldExist(field_name_counter_)){
    field_found=true;
    int int_tmp;
    ifile.scanField(field_name_counter_,int_tmp);
    counter=(unsigned int) int_tmp;
  }
  return field_found;
}


bool CounterBase::isCounterFieldInFile(IFile& ifile) {
  return ifile.FieldExist(field_name_counter_);
}


void CounterBase::writeCounterInfoToFile(OFile& ofile) const {
  if(isActive){
    // ofile.addConstantField(field_name_time_).printField(field_name_time_,getTimeStep()*getStep());
    ofile.addConstantField(field_name_counter_).printField(field_name_counter_,(int) counter);
  }
}


}
