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
#ifndef __PLUMED_variational_IndicesBase_h
#define __PLUMED_variational_IndicesBase_h

#include <vector>
#include <string>

namespace PLMD{


/// \ingroup TOOLBOX
class IndicesBase
{
public:
  // the type of 1D index
  typedef size_t index_t;
  // typedef unsigned int index_t;
private:
  unsigned int ndimensions_;
  std::vector<unsigned int> nelements_per_dimension_;
  index_t nelements_total_;
  std::vector<std::string> elements_descriptions_;
  std::vector<std::string> dimension_labels_;
public:
  IndicesBase();
  IndicesBase(const std::vector<unsigned int>&);
  ~IndicesBase() {}
  //
  std::vector<unsigned int> getNumberOfElementsPerDimension() const;
  unsigned int getNumberOfElementsPerDimension(const unsigned int) const;
  index_t getTotalNumberOfElements() const;
  unsigned int numberOfDimension() const;
  //
  index_t getIndex(const std::vector<unsigned int>&) const;
  std::vector<unsigned int> getIndices(const index_t) const;
  //
  std::string getElementDescription(const index_t) const;
  std::string getElementDescription(const std::vector<unsigned int>&) const;
  std::vector<std::string> getAllElementsDescriptions() const;
  void setElementDescription(const index_t, const std::string);
  void setElementDescription(const std::vector<unsigned int>&, const std::string);
  void setAllElementDescriptions(const std::string);
  void setAllElementDescriptions(const std::vector<std::string>&);
  //
  std::string getDimensionLabel(const unsigned int) const;
  std::vector<std::string> getAllDimensionLabels() const;
  void setDimensionLabel(const unsigned int, const std::string);
  void setAllDimensionLabels(const std::string);
  void setAllDimensionLabels(const std::vector<std::string>);
protected:
  void setupIndices(const std::vector<unsigned int>&);
};
}
#endif
