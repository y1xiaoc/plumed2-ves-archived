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
#ifndef __PLUMED_variational_CoeffsMatrix_h
#define __PLUMED_variational_CoeffsMatrix_h

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
class CoeffsMatrix:
  public CounterBase,
  public CoeffsBase
{
public:
private:
  //
  index_t size_;
  index_t nrows_;
  index_t ncolumns_;
  //
  std::vector<double> data;
  //
  bool symmetric_;
  bool diagonal_;
  //
  std::string output_fmt_; // format for output
  //
  void setupMatrix();
  //
public:
  CoeffsMatrix(
    const std::string&,
    const std::vector<std::string>&,
    const std::vector<unsigned int>&,
    const bool symmetric=true,
    const bool diagonal=false,
    const bool use_counter=false);
  CoeffsMatrix(
    const std::string&,
    std::vector<Value*>,
    std::vector<BasisFunctions*>,
    const bool symmetric=true,
    const bool diagonal=false,
    const bool use_counter=false);
  ~CoeffsMatrix(){}
  //
  index_t getSize() const;
  //
  bool isSymmetric() const;
  bool isDiagonal() const;
  //
  index_t getMatrixIndex(const index_t, const index_t) const;
  //
  // clear coeffs
  void clear();
  // get value
  double getValue(const index_t, const index_t) const;
  double getValue(const std::vector<unsigned int>&, const std::vector<unsigned int>&) const;
  // set value
  void setValue(const index_t, const index_t, const double);
  void setValue(const std::vector<unsigned int>&, const std::vector<unsigned int>&, const double);
  // add to value
  void addToValue(const index_t, const index_t, const double);
  void addToValue(const std::vector<unsigned int>&, const std::vector<unsigned int>&, const double);
  // scale all values
  void scaleAllValues(const double);
  // set all values
  void setValues(const double);
  // add to all values
  void addToValues(const double);
  // copy values for another Coeffs instance
  void setFromOtherCoeffsMatrix(CoeffsMatrix*);
  void setFromOtherCoeffsMatrix(CoeffsMatrix*,const double);
  void addFromOtherCoeffsMatrix(CoeffsMatrix*);
  void addFromOtherCoeffsMatrix(CoeffsMatrix*,const double);
  //
  double getMinValue() const;
  double getMaxValue() const;
  //
  void randomizeValuesGaussian(int);
  //
  // file input/output stuff
  void writeMatrixInfoToFile(OFile&);
  void writeHeaderToFile(OFile&);
  void writeDataDiagonalToFile(OFile&);
  void writeDataToFile(OFile&);
  void writeToFile(OFile&);
  void writeToFile(const std::string&, const bool append_file=false);
    // set output format
  void setOutputFmt(std::string ss){ output_fmt_=ss; }

};
}


#endif
