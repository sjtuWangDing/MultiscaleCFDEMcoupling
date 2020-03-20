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

Class
    momCoupleModel

SourceFiles
    momCoupleModel.C

\*---------------------------------------------------------------------------*/

#ifndef momCoupleModel_H
#define momCoupleModel_H

#include "fvCFD.H"
#include "cfdemCloud.H"
#include "dataExchangeModel.H"

namespace Foam {

class momCoupleModel {

protected:

  const dictionary& dict_;

  cfdemCloud& particleCloud_;

  // @brief Max fluid volume fraction to calculate exchange field
  scalar maxAlpha_;

public:

  // @brief Runtime type information
  TypeName("momCoupleModel");

  // @brief Declare runtime constructor selection table
  declareRunTimeSelectionTable(
    autoPtr,
    momCoupleModel,
    dictionary,
    (
      const dictionary& dict,
      cfdemCloud& sm
    ),
    (dict, sm)
  );

  // @brief Constructors
  momCoupleModel(const dictionary& dict, cfdemCloud& sm);

  // @brief Destructor
  virtual ~momCoupleModel();

  // @brief Selector
  static autoPtr<momCoupleModel> New(const dictionary& dict,
                                     cfdemCloud& sm,
                                     word momCoupleType);

  // Member Functions

  virtual void applyDebugSettings(bool) const {}

  // @brief implicit momentum source field
  virtual tmp<volScalarField> impMomSource() const;

  // @brief explicit momentum source field
  virtual tmp<volVectorField> expMomSource() const;

  virtual dimensionedVector returnIntegralSourceField() const;

  virtual void resetMomSourceField() const = 0;

  virtual void setSourceField(volVectorField &) const {}

};  // End of momCoupleModel

}  // End of namespace Foam

#endif  // momCoupleModel_H