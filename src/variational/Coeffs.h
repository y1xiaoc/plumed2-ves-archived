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
#ifndef __PLUMED_variational_Coeffs_h
#define __PLUMED_variational_Coeffs_h

#include <vector>
#include <string>
#include <cmath>

namespace PLMD{ 


class Value;
class IFile;
class OFile;

/// \ingroup TOOLBOX
class Coeffs
{
protected:
 std::vector<double> coeffs;
 std::vector<double> aux_coeffs;
 bool use_aux_coeffs_;
 bool use_counter_;
 std::string coeffs_label_;
 std::string coeffs_type_;
 std::vector<std::string> coeffs_description_;
 std::vector<std::string> dimension_labels_;
 std::vector<unsigned int> ncoeffs_per_dimension_;
 unsigned int ncoeffs_total_;
 unsigned int dimension_;
 unsigned int counter;
 std::string fmt_; // format for output 
public:
 Coeffs(const std::string & coeffs_label, const std::string coeffs_type,
        const std::vector<std::string> & dimension_labels,
        const std::vector<unsigned int> & ncoeffs_per_dimension,
        const std::vector<std::string> coeffs_description,
        const bool use_aux_coeffs=false, const bool use_counter=false);
 ///
 void Init(const std::string & coeffs_label, const std::string coeffs_type, 
        const std::vector<std::string> & dimension_labels,
        const std::vector<unsigned int> & ncoeffs_per_dimension, 
        const std::vector<std::string> coeffs_description,
        const bool use_aux_coeffs, const bool use_counter);

/// clear coeffs
 void clearMain();
 void clearAux();
 void clear();

/// get number of basis function 
 std::vector<unsigned int> getNumberOfCoeffsPerDimension() const;
/// get dimension
 unsigned int getDimension() const;
 
/// methods to handle indices 
 std::vector<unsigned int> getIndices(const unsigned int) const;
 unsigned int getIndex(const std::vector<unsigned int> &) const;

/// get size
 unsigned getSize() const;
/// get value
 double getValue(const unsigned int index) const;
 double getValue(const std::vector<unsigned int>&) const;
 double getAuxValue(const unsigned int index) const;
 double getAuxValue(const std::vector<unsigned int>&) const;
/// set value 
 void setValue(const unsigned int, const double);
 void setValue(const std::vector<unsigned int>&, const double);
 void setAuxValue(const unsigned, const double);
 void setAuxValue(const std::vector<unsigned int>&, const double);
/// add to value
 void addValue(const unsigned int , const double); 
 void addValue(const std::vector<unsigned int>&, const double);
 void addAuxValue(const unsigned int, const double); 
 void addAuxValue(const std::vector<unsigned int>&, const double);
/// scale all values 
 void scaleAllCoeffs(const double&);
 void scaleMainCoeffs(const double&);
 void scaleAuxCoeffs(const double&);


/// file output stuff
 void writeHeader(OFile&);
 void writeToFile(OFile&,const bool);
/// set output format
 void setOutputFmt(std::string ss){fmt_=ss;}

/// counter stuff 
 void resetCounter();
 void increaseCounter();
 void addToCounter(unsigned int);
 void setCounter(unsigned int);
 unsigned int getCounter() const;

/// coeffs description stuff
void setCoeffDescription(const unsigned int, const std::string);
void setCoeffDescription(const std::vector<unsigned int>&, const std::string);
std::string getCoeffDescription(const unsigned int) const;
std::string getCoeffDescription(const std::vector<unsigned int>&) const;

 ~Coeffs(){}

};
}

  
#endif
