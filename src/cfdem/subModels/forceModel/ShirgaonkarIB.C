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
  ShirgaonkarIB
\*---------------------------------------------------------------------------*/

#include "ShirgaonkarIB.H"

namespace Foam {

cfdemDefineTypeName(Archimedes)

cfdemAddToNewFunctionMap(forceModel, Archimedes)

//! \brief Constructor
ArchimedesIB::ArchimedesIB(cfdemCloud& cloud):
  forceModel(cloud),
  subPropsDict_(cloud.couplingPropertiesDict().subDict(typeName_ + "Props")),
  voidFractionFieldName_(subPropsDict_.lookupOrDefault("voidfractionFieldName", "voidfractionNext")),
  gravityFieldName_(subPropsDict_.lookupOrDefault("gravityFieldName", "g")),
  voidFraction_(cloud.mesh().lookupObject<volScalarField>(voidFractionFieldName_)),
#if defined(version21)
  g_(cloud.mesh().lookupObject<uniformDimensionedVectorField>(gravityFieldName_))
#elif defined(version16ext) || defined(version15)
  g_(dimensionedVector(cloud.mesh().lookupObject<IOdictionary>("environmentalProperties").lookup(environmentalProperties)).value())
#endif
{

}

ArchimedesIB::~ArchimedesIB() {}

void ArchimedesIB::setForce() {
  
}

} // namespace Foam
