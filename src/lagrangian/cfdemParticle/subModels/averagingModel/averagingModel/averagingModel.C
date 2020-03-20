/*---------------------------------------------------------------------------*\
    CFDEMcoupling - Open Source CFD-DEM coupling

    CFDEMcoupling is part of the CFDEMproject
    www.cfdem.com
                                Christoph Goniva, christoph.goniva@cfdem.com
                                Copyright 2009-2012 JKU Linz
                                Copyright 2012-     DCS Computing GmbH, Linz
-------------------------------------------------------------------------------
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
\*---------------------------------------------------------------------------*/

#include "averagingModel.H"

namespace Foam {

defineTypeNameAndDebug(averagingModel, 0);

defineRunTimeSelectionTable(averagingModel, dictionary);

averagingModel::averagingModel(const dictionary& dict, cfdemCloud& sm):
  dict_(dict),
  particleCloud_(sm),
  UsWeightField_
  (
    IOobject
    (
      "UsWeightField_",
      particleCloud_.mesh().time().timeName(),
      particleCloud_.mesh(),
      IOobject::NO_READ,
      IOobject::AUTO_WRITE
    ),
    particleCloud_.mesh(),
    dimensionedScalar("zero", dimensionSet(0, 0, 0, 0, 0), 0.0)
  ),
  UsPrev_
  (
    IOobject
    (
      "UsPrev",
      sm.mesh().time().timeName(),
      sm.mesh(),
      IOobject::READ_IF_PRESENT,
      IOobject::AUTO_WRITE
    ),
    sm.mesh().lookupObject<volVectorField>("Us")
    // sm.mesh(),
    // dimensionedVector("zero", dimensionSet(0, 1, -1, 0, 0), vector::zero)
  ),
  UsNext_
  (
    IOobject
    (
      "UsNext",
      sm.mesh().time().timeName(),
      sm.mesh(),
      IOobject::READ_IF_PRESENT,
      IOobject::AUTO_WRITE
    ),
    sm.mesh().lookupObject<volVectorField>("Us")
    // sm.mesh(),
    // dimensionedVector("zero", dimensionSet(0, 1, -1, 0, 0), vector::zero)
  ) {}

// @brief Destructor
averagingModel::~averagingModel() {}

void averagingModel::applyDebugSettings(bool debug) const {
  if(!debug) {
    UsWeightField_.writeOpt() = IOobject::NO_WRITE;
    UsPrev_.writeOpt() = IOobject::NO_WRITE;
    UsNext_.writeOpt() = IOobject::NO_WRITE;
  }
}

void averagingModel::setVectorSum(volVectorField& field,
                                  double**& value,
                                  double**& weight,
                                  double** const& mask) const {
  label cellI(-1);
  vector valueVec(0, 0, 0);
  for(int index = 0; index < particleCloud_.numberOfParticles(); index++) {
    for(int subCell = 0; subCell < particleCloud_.cellsPerParticle()[index][0]; subCell++) {
      cellI = particleCloud_.cellIDs()[index][subCell];
      if (cellI >= 0) {
        for (int i = 0; i < 3; ++i) { valueVec[i] = value[index][i]; }
        field[cellI] += valueVec * weight[index][subCell];
      }
    }
  }
  // correct cell values to patches
  field.correctBoundaryConditions();
}

void averagingModel::setScalarSum(volScalarField& field,
                                  double**& value,
                                  double** const& weight,
                                  double** const& mask) const {
  label cellI(-1);
  scalar valueScal(0);
  for(int index = 0; index < particleCloud_.numberOfParticles(); index++) {
    for(int subCell = 0; subCell < particleCloud_.cellsPerParticle()[index][0]; subCell++) {
      cellI = particleCloud_.cellIDs()[index][subCell];
      if (cellI >= 0) {
        valueScal = value[index][0];
        field[cellI] += valueScal * weight[index][subCell];
      }
    }
  }
  // correct cell values to patches
  field.correctBoundaryConditions();
}

tmp<volVectorField> averagingModel::UsInterp() const {
  scalar tsf = particleCloud_.dataExchangeM().timeStepFraction();
  // Info << "using Us blend, tsf=" << tsf << endl;
  return tmp<volVectorField>(new volVectorField("Us_averagingModel", (1 - tsf) * UsPrev_ + tsf * UsNext_));
}

}  // End of namespace Foam
