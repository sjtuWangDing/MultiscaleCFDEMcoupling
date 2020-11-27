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
  The IB voidFraction model is supposed to be used when a particle
  (or its representation) is bigger than a CFD cell. The voidfraction field
  is set in those cell whose centres are inside the particle. The model is
  specially designed for cfdemSolverIB and creates a smooth transition of
  the voidfraction at the particle surface. Cells which are only partially
  covered by solid are marked by voidfraction values between 0 and 1
  respectively.

  The region of influence of a particle can be increased artificially
  by “scaleUpVol”, which blows up the particles, but keeps their volume
  (for voidfraction calculation) constant.

Syntax
  voidfractionModel IB;
  IBProps
  {
    maxCellsPerParticle number1;
    alphaMin number2;
    scaleUpVol number3;
    checkPeriodicCells ;
  }
  - number1 = maximum number of cells covered by a particle (search will
              fail when more than number1 cells are covered by the particle)
  - number2 = minimum limit for voidfraction
  - number3 = diameter of the particle’s representation is artificially
              increased according to number3 * Vparticle, volume remains
              unaltered!
  - checkPeriodicCells = (optional, default false) flag for considering the
              minimal distance to all periodic images of this particle

Class
  IBVoidFraction
\*---------------------------------------------------------------------------*/

#ifndef __IB_VOID_FRACTION_H__
#define __IB_VOID_FRACTION_H__

#include "subModels/voidFractionModel/voidFractionModel.H"

namespace Foam {

class IBVoidFraction: public voidFractionModel {

public:

  cfdemTypeName("IB")

  cfdemDefineNewFunctionAdder(voidFractionModel, IBVoidFraction)

  //! \brief Constructor
  IBVoidFraction(cfdemCloud& cloud);

  //! \brief Destructor
  virtual ~IBVoidFraction();

protected:

  /*!
   * \brief 设置单个颗粒的体积分数场
   * \param index 颗粒索引
   */
  void setVolumeFractionForSingleParticle(int index) const;

  /*!
   * \brief 计算距离系数
   * \param index           <[in] 颗粒索引
   * \param pointInside     <[in] 网格中心
   * \param pointOutside    <[in] 网格角点
   */
  double segmentParticleIntersection(int index,
                                     const Foam::vector& pointInside,
                                     const Foam::vector& pointOutside) const;

};

} // namespace Foam

#endif // __IB_VOID_FRACTION_H__
