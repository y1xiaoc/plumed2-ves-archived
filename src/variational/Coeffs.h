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
class BasisFunctions;

/// \ingroup TOOLBOX
class Coeffs
{
protected:
 std::vector<double> coeffs;
 std::vector<double> aux_coeffs;
 bool useaux_;
 bool usecounter_;
 // if the coeffs are for a linear basis set expansion composed of products of 1-D basis functions
 bool isbasisfcoeffs_;
 std::string coeffs_label_;
 std::string coeffs_type_;
 std::vector<std::string> coeffs_descriptions_;
 std::vector<std::string> dimension_labels_;
 std::vector<unsigned int> ncoeffs_per_dimension_;
 unsigned int ncoeffs_total_;
 unsigned int dimension_;
 unsigned int counter;
 std::string fmt_; // format for output 
 //
 std::vector<std::string> basisf_type_;
 std::vector<unsigned int> basisf_order_;
 std::vector<unsigned int> basisf_size_;
 std::vector<double> basisf_min_;
 std::vector<double> basisf_max_;
 std::vector<std::string> basisf_keywords_;
 //
 void Init(const std::string& coeffs_label, const std::string& coeffs_type, 
        const std::vector<std::string>& dimension_labels,
        const std::vector<unsigned int>& ncoeffs_per_dimension, 
        const bool use_aux_coeffs, const bool use_counter);
 void setupBasisFunctionsInfo(std::vector<BasisFunctions*>);
 void setupBasisFunctionFromFile(const std::vector<std::string>&);
public:
 Coeffs(const std::string& coeffs_label, const std::string& coeffs_type,
        const std::vector<std::string>& dimension_labels,
        const std::vector<unsigned int>& ncoeffs_per_dimension,
        const bool use_aux_coeffs=false, const bool use_counter=false);
 //
 Coeffs(const std::string& coeffs_label,
        std::vector<Value*> args,
        std::vector<BasisFunctions*> basisf,
        const bool use_aux_coeffs=false, const bool use_counter=false);
 // clear coeffs
 void clearMain();
 void clearAux();
 void clear();

 std::string getLabel() const;
 std::string getType() const;
 bool isBasisFunctionCoeffs() const;
 bool hasAuxCoeffs() const;
 bool hasCounter() const;
 std::vector<unsigned int> getNumberOfCoeffsPerDimension() const;
 unsigned int getSize() const;
 unsigned int getDimension() const;
 
/// methods to handle indices 
 std::vector<unsigned int> getIndices(const unsigned int) const;
 unsigned int getIndex(const std::vector<unsigned int> &) const;

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
 void setValueAndAux(const unsigned, const double, const double);
 void setValueAndAux(const std::vector<unsigned int>&, const double, const double);
/// add to value
 void addToValue(const unsigned int , const double); 
 void addToValue(const std::vector<unsigned int>&, const double);
 void addToAuxValue(const unsigned int, const double); 
 void addToAuxValue(const std::vector<unsigned int>&, const double);
/// scale all values 
 void scaleAllValues(const double);
 void scaleOnlyMainValues(const double);
 void scaleOnlyAuxValues(const double);
 /// set all values 
 void setValues(const double);
 void setAuxValues(const double);
/// add to all values
 void addToValues(const double);
 void addToAuxValues(const double);
/// set Aux values equal to main
 void setAuxEqualToMain();
/// copy values for another Coeffs instance
 void copyOtherCoeffs(Coeffs*);
 void copyAndScaleOtherCoeffs(Coeffs*, const double);
/// Random coeffs
 void randomizeCoeffs();


/// file output stuff
 void writeHeader(OFile&);
 void writeToFile(OFile&,const bool print_description=false);
 void writeToFile(const std::string&,const bool print_description=false, const bool append_file=false);
 unsigned int readFromFile(IFile&, const bool ignore_missing_coeffs=false);
 unsigned int readFromFile(const std::string&, const bool ignore_missing_coeffs=false);
 static Coeffs* createFromFile(IFile&, const bool ignore_missing_coeffs=false);
 static Coeffs* createFromFile(const std::string&, const bool ignore_missing_coeffs=false);
 static std::vector<std::string> getBasisFunctionKeywordsFromFile(IFile&);
 static std::vector<std::string> getBasisFunctionKeywordsFromFile(const std::string&);
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
std::vector<std::string> getAllCoeffsDescriptions() const;
void setupCoeffsDescriptions(const std::string);

 ~Coeffs(){}

};
}

  
#endif
