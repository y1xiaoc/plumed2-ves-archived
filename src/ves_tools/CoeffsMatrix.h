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
#ifndef __PLUMED_ves_tools_CoeffsMatrix_h
#define __PLUMED_ves_tools_CoeffsMatrix_h

#include <vector>
#include <string>
#include <cmath>

#include "CoeffsBase.h"
#include "CounterBase.h"

namespace PLMD{

class Action;
class Value;
class IFile;
class OFile;
class Communicator;
class BasisFunctions;
class CoeffsVector;


/// \ingroup TOOLBOX
class CoeffsMatrix:
  public CounterBase,
  public CoeffsBase
{
public:
private:
  //
  Communicator& mycomm;
  //
  size_t size_;
  size_t nrows_;
  size_t ncolumns_;
  //
  std::vector<double> data;
  //
  bool diagonal_;
  //
  std::string output_fmt_; // format for output
  //
  void setupMatrix();
  //
public:
  explicit CoeffsMatrix(
    const std::string&,
    const std::vector<std::string>&,
    const std::vector<unsigned int>&,
    Communicator& cc,
    const bool diagonal=true,
    const bool use_counter=false);
  //
  explicit CoeffsMatrix(
    const std::string&,
    std::vector<Value*>&,
    std::vector<BasisFunctions*>&,
    Communicator& cc,
    const bool diagonal=true,
    const bool use_counter=false);
  //
  explicit CoeffsMatrix(
    const std::string&,
    std::vector<std::vector<Value*> >& argsv,
    std::vector<std::vector<BasisFunctions*> >& basisfv,
    Communicator& cc,
    const bool diagonal=true,
    const bool use_counter=false);
  //
  explicit CoeffsMatrix(
    const std::string&,
    CoeffsVector*,
    Communicator& cc,
    const bool diagonal=true,
    const bool use_counter=false);
  //
  ~CoeffsMatrix(){}
  //
  size_t getSize() const;
  //
  bool isSymmetric() const;
  bool isDiagonal() const;
  //
  bool sameShape(CoeffsVector&) const;
  bool sameShape(CoeffsMatrix&) const;
  //
  void sumCommMPI();
  void sumCommMPI(Communicator&);
  //
  void sumMultiSimCommMPI(Communicator&);
  //
  size_t getMatrixIndex(const size_t, const size_t) const;
  //
  // clear coeffs
  void clear();
  // get value
  double getValue(const size_t, const size_t) const;
  double getValue(const std::vector<unsigned int>&, const std::vector<unsigned int>&) const;
  // set value
  void setValue(const size_t, const size_t, const double);
  void setValue(const std::vector<unsigned int>&, const std::vector<unsigned int>&, const double);
  double& operator()(const size_t, const size_t);
  const double& operator()(const size_t, const size_t) const;
  double& operator()(const std::vector<unsigned int>&, const std::vector<unsigned int>&);
  const double& operator()(const std::vector<unsigned int>&, const std::vector<unsigned int>&) const;
  //
  friend CoeffsVector operator*(const CoeffsMatrix&, const CoeffsVector&);
  // add to value
  void addToValue(const size_t, const size_t, const double);
  void addToValue(const std::vector<unsigned int>&, const std::vector<unsigned int>&, const double);
  // scale all values
  void scaleAllValues(const double);
  CoeffsMatrix& operator*=(const double);
  friend CoeffsMatrix operator*(const double, const CoeffsMatrix&);
  friend CoeffsMatrix operator*(const CoeffsMatrix&, const double);
  CoeffsMatrix& operator*=(const CoeffsMatrix&);
  CoeffsMatrix operator*(const CoeffsMatrix&) const;
  // set all values
  void setValues(const double);
  void setValues(const std::vector<double>&);
  void setValues(const CoeffsMatrix&);
  CoeffsMatrix& operator=(const double);
  CoeffsMatrix& operator=(const std::vector<double>&);
  CoeffsMatrix& operator=(const CoeffsMatrix&);
  // add to all values
  CoeffsMatrix operator+() const;
  CoeffsMatrix operator-() const;
  void addToValues(const double);
  void addToValues(const std::vector<double>&);
  void addToValues(const CoeffsMatrix&);
  void subtractFromValues(const double);
  void subtractFromValues(const std::vector<double>&);
  void subtractFromValues(const CoeffsMatrix&);
  CoeffsMatrix& operator+=(const double);
  friend CoeffsMatrix operator+(const double, const CoeffsMatrix&);
  friend CoeffsMatrix operator+(const CoeffsMatrix&, const double);
  CoeffsMatrix& operator+=(const std::vector<double>&);
  friend CoeffsMatrix operator+(const std::vector<double>&, const CoeffsMatrix&);
  friend CoeffsMatrix operator+(const CoeffsMatrix&, const std::vector<double>&);
  CoeffsMatrix& operator-=(const double);
  friend CoeffsMatrix operator-(const double, const CoeffsMatrix&);
  friend CoeffsMatrix operator-(const CoeffsMatrix&, const double);
  CoeffsMatrix& operator-=(const std::vector<double>&);
  friend CoeffsMatrix operator-(const std::vector<double>&, const CoeffsMatrix&);
  friend CoeffsMatrix operator-(const CoeffsMatrix&, const std::vector<double>&);
  CoeffsMatrix& operator+=(const CoeffsMatrix&);
  CoeffsMatrix operator+(const CoeffsMatrix&) const;
  CoeffsMatrix& operator-=(const CoeffsMatrix&);
  CoeffsMatrix operator-(const CoeffsMatrix&) const;
  //
  double getMinValue() const;
  double getMaxValue() const;
  //
  void randomizeValuesGaussian(int);
  //
  // file input/output stuff
  void writeToFile(OFile&, const double current_time=-1.0);
  void writeToFile(const std::string&, const double current_time=-1.0, const bool append_file=false, Action* action_ptr=NULL);
private:
  void writeDataToFile(OFile&);
  void writeMatrixInfoToFile(OFile&);
  void writeHeaderToFile(OFile&, const double current_time=-1.0);
  void writeDataDiagonalToFile(OFile&);
  void writeDataFullToFile(OFile&);
public:
    // set output format
  void setOutputFmt(std::string ss){ output_fmt_=ss; }
  std::string getOutputFmt() const {return output_fmt_;}
  Communicator& getCommunicator() const {return mycomm;}

};
}


#endif
