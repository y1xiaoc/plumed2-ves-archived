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
#ifndef __PLUMED_variational_CoeffsVector_h
#define __PLUMED_variational_CoeffsVector_h

#include <vector>
#include <string>
#include <cmath>

#include "CoeffsBase.h"
#include "CounterBase.h"

namespace PLMD{


class Value;
class IFile;
class OFile;
class BasisFunctions;

/// \ingroup TOOLBOX
class CoeffsVector:
  public CounterBase,
  public CoeffsBase
{
public:
private:
  std::vector<double> data;
  std::string output_fmt_; // format for output
  //
public:
  CoeffsVector(
    const std::string&,
    const std::vector<std::string>&,
    const std::vector<unsigned int>&,
    const bool use_counter=false);
  CoeffsVector(
    const std::string&,
    std::vector<Value*>,
    std::vector<BasisFunctions*>,
    const bool use_counter=false);
  ~CoeffsVector(){}
  //
  index_t getSize() const;
    // clear coeffs
  void clear();
  //
  bool sameShape(const CoeffsVector other_coeffsvector) const;
  // get value
  double getValue(const index_t) const;
  double getValue(const std::vector<unsigned int>&) const;
  double& operator[](const index_t index);
  const double& operator[](const index_t index) const;
  // set value
  void setValue(const index_t, const double);
  void setValue(const std::vector<unsigned int>&, const double);
    // add to value
  void addToValue(const index_t, const double);
  void addToValue(const std::vector<unsigned int>&, const double);
  // scale all values
  void scaleAllValues(const double);
  CoeffsVector operator*=(const double);
  friend CoeffsVector operator*(const double, const CoeffsVector&);
  friend CoeffsVector operator*(const CoeffsVector&, const double);
  CoeffsVector operator*=(const CoeffsVector&);
  CoeffsVector operator*(const CoeffsVector&) const;
  // set all values
  void setValues(const double);
  void setValues(const std::vector<double>&);
  void setValues(const CoeffsVector&);
  CoeffsVector operator=(const double);
  CoeffsVector operator=(const std::vector<double>&);
  CoeffsVector operator=(const CoeffsVector&);
  // add to all values
  void addToValues(const double);
  void addToValues(const CoeffsVector&);
  CoeffsVector operator+=(const double);
  friend CoeffsVector operator+(const double, const CoeffsVector&);
  friend CoeffsVector operator+(const CoeffsVector&, const double);
  CoeffsVector operator+=(const CoeffsVector&);
  CoeffsVector operator+(const CoeffsVector&) const;

  //
  double getMinValue() const;
  double getMaxValue() const;
  double getNorm() const;
  //
  void normalizeCoeffs();
  // Random values
  void randomizeValuesGaussian(int);

  // file input/output stuff
  void writeToFile(const std::string&, const bool print_description=false, const bool append_file=false);
  void writeToFile(OFile&,const bool print_description=false);
  static void writeToFile(const std::string&, const std::vector<CoeffsVector>&, const bool print_description=false, const bool append_file=false);
  static void writeToFile(OFile&, const std::vector<CoeffsVector>&, const bool print_description=false);
private:
  void writeHeaderToFile(OFile&) const;
  static void writeDataToFile(OFile&, const std::vector<CoeffsVector>&, const bool print_description=false);
public:
  unsigned int readFromFile(IFile&, const bool ignore_missing_coeffs=false, const bool ignore_coeffs_info=false);
  unsigned int readFromFile(const std::string&, const bool ignore_missing_coeffs=false, const bool ignore_coeffs_info=false);
private:
  void readHeaderFromFile(IFile&, const bool ignore_coeffs_info=false);
  unsigned int readDataFromFile(IFile&, const bool ignore_missing_coeffs=false);
public:
  // set output format
  void setOutputFmt(std::string ss){ output_fmt_=ss; }
  std::string getOutputFmt() const { return output_fmt_; }

};
}


#endif