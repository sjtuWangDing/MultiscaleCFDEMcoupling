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
  demoModel
\*---------------------------------------------------------------------------*/

#include "voidFractionModel.H"

namespace Foam {

cfdemDefineTypeName(voidFractionModel)

cfdemDefineNewFunctionMap(voidFractionModel)

cfdemDefineConstructNewFunctionMap(voidFractionModel)

cfdemDefineDestroyNewFunctionMap(voidFractionModel)

cfdmeDefineBaseTypeNew(autoPtr, voidFractionModel, (cfdemCloud& cloud, const dictionary& dict), dict, (cloud))

//! \brief Constructor
voidFractionModel::voidFractionModel(cfdemCloud& cloud)
  : cloud_(cloud),
    voidFractionPrev_(
      IOobject(
        "voidFractionPrev",
        cloud.mesh().time().timeName(),
        cloud.mesh(),
        IOobject::READ_IF_PRESENT, // or MUST_READ,
        IOobject::AUTO_WRITE
      ),
      // cloud.mesh().lookupObject<volScalarField>("voidFraction")
      cloud.mesh(),
      dimensionedScalar("zero", dimensionSet(0, 0, 0, 0, 0), 1)
    ),
    voidFractionNext_(
      IOobject(
        "voidFractionNext",
        cloud.mesh().time().timeName(),
        cloud.mesh(),
        IOobject::READ_IF_PRESENT, // or MUST_READ,
        IOobject::AUTO_WRITE
      ),
      // cloud.mesh().lookupObject<volScalarField>("voidFraction")
      cloud.mesh(),
      dimensionedScalar("zero", dimensionSet(0, 0, 0, 0, 0), 1)
    ),
    volumeFractionPrev_(
      IOobject(
        "volumeFractionPrev",
        cloud.mesh().time().timeName(),
        cloud.mesh(),
        IOobject::READ_IF_PRESENT, // or MUST_READ,
        IOobject::AUTO_WRITE
      ),
      // cloud.mesh().lookupObject<volScalarField>("volumeFraction")
      cloud.mesh(),
      dimensionedScalar("zero", dimensionSet(0, 0, 0, 0, 0), 1)
    ),
    volumeFractionNext_(
      IOobject(
        "volumeFractionNext",
        cloud.mesh().time().timeName(),
        cloud.mesh(),
        IOobject::READ_IF_PRESENT, // or MUST_READ,
        IOobject::AUTO_WRITE
      ),
      // cloud.mesh().lookupObject<volScalarField>("volumeFraction")
      cloud.mesh(),
      dimensionedScalar("zero", dimensionSet(0, 0, 0, 0, 0), 1)
    ),
    maxCellsNumPerFineParticle_(1),
    maxCellsNumPerMiddleParticle_(1),
    maxCellsNumPerCoarseParticle_(1) {}

//! \brief Destructor
voidFractionModel::~voidFractionModel() {}



} // namespace Foam
