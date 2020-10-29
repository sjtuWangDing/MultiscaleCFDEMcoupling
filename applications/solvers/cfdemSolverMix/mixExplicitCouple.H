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
    mixExplicitCouple

SourceFiles
    mixExplicitCouple.C

\*---------------------------------------------------------------------------*/

#ifndef mixExplicitCouple_H
#define mixExplicitCouple_H

#include "momCoupleModel.H"

namespace Foam {

class mixExplicitCouple: public momCoupleModel {

private:

  dictionary propsDict_;

  //! @brief Limit for for exchange field
  vector fLimit_;

  mutable volVectorField fPrev_;

  mutable volVectorField fNext_;

  //! @brief 单位体积上颗粒对流体的作用力源项场(单位: N / m^3)
  mutable volVectorField sourceField_;

private:

  vector arrayToField(label cellI) const;

public:

  //! @brief Runtime type information
  TypeName("mixExplicitCouple");

  //! @brief Constructors
  mixExplicitCouple(const dictionary& dict, cfdemCloud& sm);

  //! @brief Destructor
  ~mixExplicitCouple();

  // Member Functions

  void applyDebugSettings(bool debug) const {
    if (!debug) {
      fPrev_.writeOpt() = IOobject::NO_WRITE;
      fNext_.writeOpt() = IOobject::NO_WRITE;
    }
  }

  void resetMomSourceField() const {
    fPrev_ == fNext_;
    fNext_ == dimensionedVector("zero", fNext_.dimensions(), vector::zero);
  }

  void setSourceField(volVectorField& field) const {
    sourceField_ = field;
  }

  //! @brief 累加当前处理器上的全部源项
  // @note 返回值单位为 N，而不是 N / m^3
  dimensionedVector returnIntegralSourceField() const;

  tmp<volVectorField> expMomSource() const;

};  // End of class mixImplicitCouple

}  // End of namespace Foam

#endif  // mixExplicitCouple_H
