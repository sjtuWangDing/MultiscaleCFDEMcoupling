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
  forceSubModel
\*---------------------------------------------------------------------------*/

#ifndef __FORCE_SUB_MODEL_H__
#define __FORCE_SUB_MODEL_H__

#include "cloud/cfdemCloud.H"

namespace Foam {

/*!
  * \brief forceSwitch
  *
  * treatForceExplicit
  *    true - 在CFD的动量方程中，耦合力都为显式力
  *    false - 在CFD的动量方程中，耦合力都为半隐式力 (default)
  *
  * treatForceDEM
  *    true - 仅在 DEM 中考虑耦合力，即仅考虑流体对颗粒的作用，不考虑颗粒对流体的反作用
  *           在 resolved 方法中，比如 cfdemSolverIB，ArchimedesIB force 会直接设置 treatForceDEM 为 true
  *           因为颗粒对流体的作用力是通过虚拟域方法得到
  *    false - 即 CFD 和 DEM 都考虑耦合力 (default)
  *
  * implForceDEM
  *    true - 在每个耦合时间步，流体的速度和阻力系数都被传递到 DEM 中，从而在每个 DEM 时间步中，
  *           使用上一个耦合时间步中的阻力系数和流体速度，与当前颗粒速度一起计算颗粒受到的阻力
  *    false - 在每个耦合时间步中，流体对颗粒的阻力被传递到 DEM 中，并且在接下来
  *            的 DEM 时间步中，这个力保持不变，直到下个耦合时间步 (default)
  *
  * verbose
  *    true - 调试信息输出到屏幕
  *    false - 调试信息不输出到屏幕 (default)
  */
enum ESwitch {
  treatForceExplicit = 0,
  treatForceDEM,
  implForceDEM,
  verbose
};

/*!
 * \brief The force sub model is designed to hold the settings a force model can have.
 *        For now it handles the treatForceExplicit, treatForceDEM and implForceDEM option.
 */
class forceSubModel {
public:
  /* --------------------------------- base class and functions -------------------------- */
  class Switches {
  public:
    static constexpr int kNum = 4;
    static constexpr const char* kNameList[] = {
      "treatForceExplicit",
      "treatForceDEM",
      "implForceDEM",
      "verbose"
    };
    Switches(): value_(0) {}
    void setTrue(ESwitch value) {
      switch(value) {
        case treatForceExplicit:
          value_ |= (1 << treatForceExplicit);
          break;
        case treatForceDEM:
          value_ |= (1 << treatForceDEM);
          break;
        case implForceDEM:
          value_ |= (1 << implForceDEM);
          break;
        case verbose:
          value_ |= (1 << verbose);
          break;
        default:
          FatalError << "Error: illegal switch enum: " << value << " in forceSubModel"
            << abort(FatalError);
      }
    }
    bool isTrue(ESwitch value) const {
      switch(value) {
        case treatForceExplicit:
          return value_ & (1 << treatForceExplicit);
        case treatForceDEM:
          return value_ & (1 << treatForceDEM);
        case implForceDEM:
          return value_ & (1 << implForceDEM);
        case verbose:
          return value_ & (1 << verbose);
        default:
          FatalError << "Error: illegal switch enum: " << value << " in forceSubModel"
            << abort(FatalError);
      }
      return false;
    }
  private:
    unsigned int value_;
  };

  //! \brief Constructor
  forceSubModel(cfdemCloud& cloud, forceModel& forceModel);

  //! \brief Destructor
  ~forceSubModel();

  /* --------------------------------- virtual functions --------------------------------- */

  /* --------------------------------- member functions ---------------------------------- */

  /*!
   * \param index                  <[in] 颗粒索引
   * \param dragTot                <[in] 索引为 index 的颗粒受到的总阻力
   * \param dragEx                 <[in] 索引为 index 的颗粒受到的显式阻力
   * \param Ufluid = vector::zero  <[in] 索引为 index 的颗粒中心处流体速度(可以指定是否使用插值模型计算)
   * \param scalar Cd = 0          <[in] 颗粒阻力系数
   */
  void partToArray(const label& index,
                   const Foam::vector& dragTot,
                   const Foam::vector& dragEx,
                   const Foam::vector& Ufluid = Foam::vector::zero,
                   scalar Cd = scalar(0)) const;

  void readSwitches(const dictionary& subPropsDict);

  /* --------------------------------- public interfaces --------------------------------- */

protected:

  cfdemCloud& cloud_;

  forceModel& forceModel_;

  Switches switches_;

  //! \brief 用户自定义动力粘度，该粘度仅用于计算颗粒所受阻力
  volScalarField nu_;
};

} // namespace Foam

#endif // __FORCE_SUB_MODEL_H__
