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
    mixImplicitCouple

SourceFiles
    mixImplicitCouple.C

\*---------------------------------------------------------------------------*/

#ifndef mixImplicitCouple_H
#define mixImplicitCouple_H

#include "momCoupleModel.H"

namespace Foam {

class mixImplicitCouple: public momCoupleModel {

private:

  dictionary propsDict_;

  //! @brief 流体速度场的名称
  word velFieldName_;

  //! @brief 局部平均颗粒速度场的名称
  word granVelFieldName_;

  //! @brief 空隙率场的名称
  word voidfractionFieldName_;

  //! @brief 速度场的常引用
  const volVectorField& U_;

  //! @brief 局部平均颗粒速度场的常引用
  const volVectorField& Us_;

  //! @brief 空隙率场的常引用
  const volScalarField& alpha_;

  //! @brief 动量交换场的最大值
  scalar KslLimit_;

  //! @brief 动量交换场
  //! @note Ksl = 单位体积上颗粒对流体的作用力 / (Us - U)，单位为 (N/m^3) / (m/s) = kg / (m^3 * s)
  mutable volScalarField KslPrev_;
  mutable volScalarField KslNext_;

public:

  //! @brief Runtime type information
  TypeName("mixImplicitCouple");

  //! @brief Constructors
  mixImplicitCouple(const dictionary& dict, cfdemCloud& sm);

  //! @brief Destructor
  ~mixImplicitCouple();

  // Member Functions

  void applyDebugSettings(bool debug) const {
    if (!debug) {
      KslPrev_.writeOpt() = IOobject::NO_WRITE;
      KslNext_.writeOpt() = IOobject::NO_WRITE;
    }
  }

  void resetMomSourceField() const{
    KslPrev_ == KslNext_;
    KslNext_ == dimensionedScalar("zero", KslNext_.dimensions(), 0.0);
  }

  tmp<volScalarField> impMomSource() const;

};  // class mixImplicitCouple

}  // namespace Foam

#endif  // mixImplicitCouple_H
