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
  This code is designed to realize coupled CFD-DEM simulations using LIGGGHTS
  and OpenFOAM(R). Note: this code is not part of OpenFOAM(R) (see DISCLAIMER).

Class
  ArchimedesIB
\*---------------------------------------------------------------------------*/

#ifndef __ARCHIMEDES_IB_H__
#define __ARCHIMEDES_IB_H__

#include "forceModel.H"

namespace Foam {

class ArchimedesIB: public forceModel {

public:

  //! @brief Runtime type information
  cfdemTypeName("ArchimedesIB")

  //! @brief Constructor
  ArchimedesIB(cfdemCloud& cloud);

  //! @brief Destructor
  ~ArchimedesIB();

  void setForce();

private:

  //! @note subPropsDict_ should be declared in front of other members
  dictionary subPropsDict_;

  std::string voidFractionFieldName_;

  std::string gravityFieldName_;

  const volScalarField& voidFraction_;

#ifdef version21
  const uniformDimensionedVectorField& g_;    // ref to gravity
#elif defined(version16ext) || defined(version15)
  const dimensionedVector& g_;    // ref to gravity
#endif
};

} // namespace Foam

#endif // __ARCHIMEDES_IB_H__
