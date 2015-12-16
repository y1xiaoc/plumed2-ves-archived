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
#ifndef __PLUMED_ves_tools_CounterBase_h
#define __PLUMED_ves_tools_CounterBase_h

#include <string>

namespace PLMD{

class IFile;
class OFile;

/// \ingroup TOOLBOX
class CounterBase{
private:
  unsigned int counter;
  std::string field_name_counter_;
  std::string field_name_time_;
  bool isActive;
public:
  explicit CounterBase(const bool active=true);
  ~CounterBase() {}
  //
  void turnOnCounter();
  void turnOffCounter();
  bool isCounterActive();
  //
  void resetCounter();
  void increaseCounter();
  void addToCounter(const unsigned int);
  void setCounter(const unsigned int);
  unsigned int getCounter() const;
  double getCounterDbl() const;
  //
  void setCounterFieldName(const std::string&);
  std::string getCounterFieldName() const;
  void setTimeFieldName(const std::string&);
  std::string getTimeFieldName() const;
  //
  bool getCounterFieldFromFile(IFile&);
  bool isCounterFieldInFile(IFile&);
  void writeCounterInfoToFile(OFile&) const;
  //
};
}

#endif
