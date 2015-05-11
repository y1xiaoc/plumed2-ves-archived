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
#include <vector>
#include <cmath>
#include <iostream>
#include <sstream>
#include <cstdio>
#include <cfloat>

#include "Coeffs.h"
#include "tools/Tools.h"
#include "core/Value.h"
#include "tools/File.h"
#include "tools/Exception.h"

using namespace std;
namespace PLMD{

Coeffs::Coeffs(const std::string & coeffs_label, const std::string coeffs_type,
        const std::vector<std::string> & dimension_labels,
        const std::vector<unsigned int> & ncoeffs_per_dimension,
        const std::vector<std::string> coeffs_description,
        const bool use_aux_coeffs, const bool use_counter)
{
 Init(coeffs_label, coeffs_type, dimension_labels, ncoeffs_per_dimension, coeffs_description, use_aux_coeffs, use_counter);
}

void Coeffs::Init(const std::string & coeffs_label, const std::string coeffs_type,
        const std::vector<std::string> & dimension_labels,
        const std::vector<unsigned int> & ncoeffs_per_dimension,
        const std::vector<std::string> coeffs_description,
        const bool use_aux_coeffs, const bool use_counter)
{
 fmt_="%14.9f"; 
 plumed_massert(ncoeffs_per_dimension.size()==dimension_labels.size(),"Coeffs: dimensions of vectors in Init(...) don't match");
 dimension_=ncoeffs_per_dimension.size(); 
 dimension_labels_=dimension_labels;
 ncoeffs_per_dimension_=ncoeffs_per_dimension;
 coeffs_label_=coeffs_label;
 coeffs_type_=coeffs_type;
 use_counter_=use_counter;
 use_aux_coeffs_=use_aux_coeffs;
 ncoeffs_total_=1;
 for(unsigned int i=0;i<dimension_;i++){ncoeffs_total_*=ncoeffs_per_dimension_[i];}
 // coeffs_description
 plumed_massert(coeffs_description.size()==ncoeffs_total_,"Coeffs: size of coeffs_descriptions don't match total number of coeffs");
 coeffs_description_=coeffs_description;
 if(use_counter_){resetCounter();}
 clear();
}

void Coeffs::clearMain(){
 coeffs.resize(ncoeffs_total_);
 for(unsigned int i=0;i<ncoeffs_total_;i++){coeffs[i]=0.0;}
}

void Coeffs::clearAux(){
 aux_coeffs.resize(ncoeffs_total_);
 for(unsigned int i=0;i<ncoeffs_total_;i++){aux_coeffs[i]=0.0;}
}

void Coeffs::clear(){
 clearMain();
 if(use_aux_coeffs_){clearAux();}
}

vector<unsigned int> Coeffs::getNumberOfCoeffsPerDimension() const {
 return ncoeffs_per_dimension_;
}

unsigned Coeffs::getSize() const {
 return ncoeffs_total_;
}

unsigned Coeffs::getDimension() const {
 return dimension_;
}

// we are flattening arrays using a column-major order
unsigned int Coeffs::getIndex(const vector<unsigned int>& indices) const {
 plumed_dbg_assert(indices.size()==dimension_);
 for(unsigned int i=0;i<dimension_;i++)
  if(indices[i]>=ncoeffs_per_dimension_[i]) {
    std::string is;
    Tools::convert(i,is);
    std::string msg="ERROR: the system is looking for a value outside the indices along the " + is;
    plumed_merror(msg+" index!");
  }
 unsigned int index=indices[dimension_-1];
 for(unsigned int i=dimension_-1;i>0;--i){
  index=index*ncoeffs_per_dimension_[i-1]+indices[i-1];
 }
 return index;
}

// we are flattening arrays using a column-major order
vector<unsigned int> Coeffs::getIndices(const unsigned int index) const {
 vector<unsigned int> indices(dimension_);
 unsigned int kk=index;
 indices[0]=(index%ncoeffs_per_dimension_[0]);
 for(unsigned int i=1;i<dimension_-1;++i){
  kk=(kk-indices[i-1])/ncoeffs_per_dimension_[i-1];
  indices[i]=(kk%ncoeffs_per_dimension_[i]);
 }
 if(dimension_>=2){
  indices[dimension_-1]=((kk-indices[dimension_-2])/ncoeffs_per_dimension_[dimension_-2]);
 }
 return indices;
}

double Coeffs::getValue(const unsigned int index) const {
 plumed_dbg_assert(index<ncoeffs_total_);
 return coeffs[index];
}

double Coeffs::getValue(const vector<unsigned int>& indices) const {
 return getValue(getIndex(indices));
}

void Coeffs::setValue(const unsigned int index, const double value){
 plumed_dbg_assert(index<ncoeffs_total_);
 coeffs[index]=value;
}

void Coeffs::setValue(const vector<unsigned int>& indices, const double value){
 setValue(getIndex(indices),value); 
}

void Coeffs::addValue(const unsigned int index, const double value){
 plumed_dbg_assert(index<ncoeffs_total_);
 coeffs[index]+=value;
}

void Coeffs::addValue(const vector<unsigned int>& indices, const double value){
 addValue(getIndex(indices),value);
}

void Coeffs::scaleAllCoeffs(const double& scalef ){
  for(unsigned int i=0;i<coeffs.size();i++){coeffs[i]*=scalef;}
  if(use_aux_coeffs_){for(unsigned int i=0;i<coeffs.size();i++){aux_coeffs[i]*=scalef;}}
}

void Coeffs::scaleMainCoeffs(const double& scalef ){
  for(unsigned int i=0;i<coeffs.size();i++){coeffs[i]*=scalef;}
}

void Coeffs::scaleAuxCoeffs(const double& scalef ){
  for(unsigned int i=0;i<aux_coeffs.size();i++){aux_coeffs[i]*=scalef;}
}

void Coeffs::writeHeader(OFile& ofile){
 ofile.addConstantField("label");
 ofile.printField("label",coeffs_label_);
 ofile.addConstantField("type");
 ofile.printField("type",coeffs_type_);
 ofile.addConstantField("ncoeffs_total");
 ofile.printField("ncoeffs_total",int(ncoeffs_total_));
 for(unsigned i=0;i<dimension_;++i){
   ofile.addConstantField("ncoeffs_" + dimension_labels_[i]);
   ofile.printField("ncoeffs_" + dimension_labels_[i],int(ncoeffs_per_dimension_[i]));
 }
 if(use_counter_){ofile.addConstantField("iteration").printField("iteration",int(counter));}
}

void Coeffs::writeToFile(OFile& ofile, const bool print_description=false){
 std::string int_fmt="%8d";
 std::string str_seperate="#!-------------------";
 std::string ilabels_prefix="idx_";
 char* s1 = new char[20];
 std::vector<unsigned int> indices;
 std::vector<std::string> ilabels(dimension_);
 for(unsigned int k=0; k<dimension_; k++){ilabels[k]=ilabels_prefix+dimension_labels_[k];}

 writeHeader(ofile); 
 for(unsigned i=0;i<ncoeffs_total_;i++){
  indices=getIndices(i);
  for(unsigned int k=0; k<dimension_; k++)
  {
   sprintf(s1,int_fmt.c_str(),indices[k]); ofile.printField(dimension_labels_[k],s1);
  }
  ofile.fmtField(" "+fmt_).printField("coeff",coeffs[i]);
  if(use_aux_coeffs_){ofile.fmtField(" "+fmt_).printField("aux_coeff",aux_coeffs[i]);}
  if(print_description){ofile.printField("description",coeffs_description_[i]);}
  sprintf(s1,int_fmt.c_str(),i); ofile.printField("index",s1);

  ofile.printField();
 }
 ofile.fmtField();
 // blank line between iterations to allow proper plotting with gnuplot
 ofile.printf("%s\n",str_seperate.c_str());
 ofile.printf("\n");
 ofile.printf("\n");
 delete [] s1;
}

// counter stuff
void Coeffs::resetCounter(){counter=0;}
void Coeffs::increaseCounter(){counter=+1;}
void Coeffs::addToCounter(unsigned int value){counter=+value;}
void Coeffs::setCounter(unsigned int value){counter=value;}
unsigned int Coeffs::getCounter() const {return counter;}

// coeffs description stuff
void Coeffs::setCoeffDescription(const unsigned int index, const std::string description)
{
 coeffs_description_[index]=description;
}

void Coeffs::setCoeffDescription(const std::vector<unsigned int>& indices, const std::string description)
{
 setCoeffDescription(getIndex(indices), description);
}

std::string Coeffs::getCoeffDescription(const unsigned int index) const 
{
 return coeffs_description_[index];
}

std::string Coeffs::getCoeffDescription(const std::vector<unsigned int>& indices) const 
{
 return getCoeffDescription(getIndex(indices));
}




}
