/*---------------------------------------------------------------------------*\
  CFDEMcoupling - Open Source CFD-DEM coupling

  CFDEMcoupling is part of the CFDEMproject
  www.cfdem.com
                              Christoph Goniva, christoph.goniva@cfdem.com
                              Copyright 2009-2012 JKU Linz
                              Copyright 2012-     DCS Computing GmbH, Linz
------------------------------------------------------------------------------
License
  This file is part of CFDEMcoupling.

  CFDEMcoupling is free software; you can redistribute it and/or modify it
  under the terms of the GNU General Public License as published by the
  Free Software Foundation; either version 3 of the License, or (at your
  option) any later version.

  CFDEMcoupling is distributed in the hope that it will be useful, but WITHOUT
  ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
  FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
  for more details.

  You should have received a copy of the GNU General Public License
  along with CFDEMcoupling; if not, write to the Free Software Foundation,
  Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

Description
  The locateModel “engine” locates the CFD cell and cellID corresponding
  to a given position. 

Syntax
  locateModel engine;
  engineProps
  {
    treeSearch switch1;
  }

Class
  engineSearch
\*---------------------------------------------------------------------------*/

#include "subModels/locateModel/engineSearch.H"

namespace Foam {

cfdemDefineTypeName(engineSearch)

cfdemAddToNewFunctionMap(locateModel, engineSearch)

//! \brief Constructor
engineSearch::engineSearch(cfdemCloud& cloud, const std::string& derivedTypeName):
  locateModel(cloud),
  subPropsDict_(
    cloud.couplingPropertiesDict().subDict((derivedTypeName.empty() ? typeName_ : derivedTypeName) + "Props")),
  treeSearch_(subPropsDict_.lookupOrDefault<bool>("treeSearch", true)),
#if defined(version30)
  searchEngine_(cloud.mesh(), polyMesh::FACE_PLANES)
#elif defined(version21)
  searchEngine_(cloud.mesh(), polyMesh::FACEPLANES)
#elif defined(version16ext)
  searchEngine_(cloud.mesh(), false)
#endif
  {}

  //! \brief Destructor
engineSearch::~engineSearch() {}

/*!
 * \brief use search engine to get cell id of particle center
 * \param numberOfParticles 颗粒数量
 * \param cellIDs   <[in] 颗粒覆盖网格的编号
 */
void engineSearch::findCell(int numberOfParticles,
                            const ::base::CITensor2& cellIDs) const {
  CHECK_EQ(numberOfParticles, static_cast<int>(cellIDs.size(0)))
    << "Number of particle is not match";
  for (index_t i = 0; i < numberOfParticles; ++i) {
    if (cloud_.getRadius(i) > Foam::SMALL) {
      Foam::vector particleCenterPos = cloud_.getPosition(i);
      // particleCenterPos: location vector
      // -1: old cell id
      // treeSearch_: use octree search
      cellIDs[i][0] = searchEngine_.findCell(particleCenterPos, -1, treeSearch_);
    } else {
      cellIDs[i][0] = -1;
    }
  }
}

label engineSearch::findSingleCell(const Foam::vector& position, label oldCellID) const {
  return searchEngine_.findCell(position, oldCellID, treeSearch_);
}

} // namespace Foam
