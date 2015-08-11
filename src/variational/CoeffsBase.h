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
#ifndef __PLUMED_variational_CoeffsBase_h
#define __PLUMED_variational_CoeffsBase_h

#include <vector>
#include <string>

namespace PLMD{


/// \ingroup TOOLBOX
class CoeffsBase
{
public:
  // the type of 1D index for the coeffs
  typedef size_t index_t;
  // typedef unsigned int index_t;
protected:
  unsigned int dimension_;
  std::vector<unsigned int> ncoeffs_per_dimension_;
  index_t ncoeffs_total_;
  std::vector<std::string> coeffs_descriptions_;
  std::vector<std::string> dimension_labels_;
public:
  CoeffsBase();
  CoeffsBase(const std::vector<unsigned int>&);
  ~CoeffsBase(){}
  //
  std::vector<unsigned int> getNumberOfCoeffsPerDimension() const;
  index_t getNumberOfCoeffs() const;
  unsigned int getDimension() const;
  //
  index_t getIndex(const std::vector<unsigned int>&) const;
  std::vector<unsigned int> getIndices(const index_t) const;
  //
  void setCoeffDescription(const index_t, const std::string);
  void setCoeffDescription(const std::vector<unsigned int>&, const std::string);
  std::string getCoeffDescription(const index_t) const;
  std::string getCoeffDescription(const std::vector<unsigned int>&) const;
  std::vector<std::string> getAllCoeffsDescriptions() const;
  void setAllCoeffDescriptions(const std::string);
  void setAllCoeffDescriptions(const std::vector<std::string>&);
protected:
  void setupIndices(const std::vector<unsigned int>&);
};
}
#endif
