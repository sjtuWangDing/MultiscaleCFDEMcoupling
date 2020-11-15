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
  voidFractionModel
\*---------------------------------------------------------------------------*/

#ifndef __VOID_FRACTION_MODEL_H__
#define __VOID_FRACTION_MODEL_H__

#include "cloud/cfdemCloud.H"
#include "base/runTimeSelectionTables.H"

namespace Foam {

//! \brief 空隙率模型
class voidFractionModel {

public:
  /* --------------------------------- base class and functions -------------------------- */

  //! \brief Runtime type information
  cfdemTypeName("voidFractionModel")

  //! \brief Declare runtime constructor selection
  cfdemDeclareRunTimeSelection(autoPtr, voidFractionModel, (cfdemCloud& cloud), (cloud))

  //! \brief Selector
  static autoPtr<voidFractionModel> New(cfdemCloud& cloud, const dictionary& dict);

  //! \brief Constructor
  voidFractionModel(cfdemCloud& cloud);

  //! \brief Destructor
  virtual ~voidFractionModel();

  /* --------------------------------- virtual functions --------------------------------- */

  /* --------------------------------- member functions ---------------------------------- */

  //! \brief 重置空隙率
  inline void resetVoidFraction() {
    voidFractionPrev_ == voidFractionNext_;
    voidFractionNext_ == dimensionedScalar("one", voidFractionNext_.dimensions(), 1.);
  }

  //! \brief 重置体积分数
  inline void resetVolumeFraction() {
    volumeFractionPrev_ == volumeFractionNext_;
    volumeFractionNext_ == dimensionedScalar("one", volumeFractionNext_.dimensions(), 1.);
  }

  /* --------------------------------- public interfaces --------------------------------- */

protected:

  cfdemCloud& cloud_;

  //! \brief 小颗粒空隙率
  volScalarField voidFractionPrev_;

  volScalarField voidFractionNext_;

  //! \brief 大颗粒体积分数
  mutable volScalarField volumeFractionPrev_;

  mutable volScalarField volumeFractionNext_;
};

} // namespace Foam

#endif // __VOID_FRACTION_MODEL_H__
