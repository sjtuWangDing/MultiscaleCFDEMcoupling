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
  cfdemCloudIB derived from cfdemCloud

Class
  Foam::cfdemCloudIB
\*---------------------------------------------------------------------------*/

#include <mutex> // std::call_once
#include "cloud/cfdemCloudIB.H"
#include "subModels/dataExchangeModel/dataExchangeModel.H"

namespace Foam {

//! \brief Constructed from mesh
cfdemCloudIB::cfdemCloudIB(const fvMesh& mesh)
  : cfdemCloud(mesh),
    meshHasUpdated_(false) {
}

//! \brief Destructor
cfdemCloudIB::~cfdemCloudIB() {}

/*!
 * \brief 更新函数
 * \note used for cfdemSolverIB
 * \param volumeFraction  <[in, out] 大颗粒体积分数
 * \param interface       <[in, out] 界面场，用于 dynamic mesh
 */
void cfdemCloudIB::evolve(volScalarField& volumeFraction,
                          volScalarField& interface) {
  Info << __func__ << ", used for cfdemSolverIB..." << endl;
  // 检查当前流体时间步是否同时也是耦合时间步
  if (dataExchangeM().checkValidCouplingStep()) {
    // 创建用于记录 coupling time step counter
    auto pCounter = std::make_shared<dataExchangeModel::CouplingStepCounter>(dataExchangeM());
    // couple(): run liggghts command and get number of particle
    setNumberOfParticles(dataExchangeM().couple());
    // realloc memory
    reAlloc();
    // 获取 DEM 数据
    getDEMData();
    for (int i = 0; i < numberOfParticles(); ++i) {
      Pout << "radius: " << radii()[i] << endl;
      Pout << "position: " << positions()[i][0] << ", " << positions()[i][1] << ", " << positions()[i][2] << endl;
      Pout << "velocity: " << velocities()[i][0] << ", " << velocities()[i][1] << ", " << velocities()[i][2] << endl;
    }
  }
  Info << __func__ << " - done\n" << endl;
}

void cfdemCloudIB::reAlloc() {
  if (numberOfParticlesChanged()) {
    int number = numberOfParticles();
    dataExchangeM().reAlloc(pCloud_.radii(), pCloud_.pRadii(), number, 1);
    dataExchangeM().reAlloc(pCloud_.positions(), pCloud_.pPositions(), number, 3);
    dataExchangeM().reAlloc(pCloud_.velocities(), pCloud_.pVelocities(), number, 3);
  }
}

void cfdemCloudIB::getDEMData() {
  dataExchangeM().getData("radius", "scalar-atom", pCloud_.pRadii());
  dataExchangeM().getData("x", "vector-atom", pCloud_.pPositions());
  dataExchangeM().getData("v", "vector-atom", pCloud_.pVelocities());
}

//! @brief 确定颗粒周围细化网格的区域(每个方向的尺寸都是颗粒尺寸的两倍)
void cfdemCloudIB::setInterface(volScalarField& interface,
                                volScalarField& refineMeshKeepStep) const {
  // 确保 call_once 函数中的 lambda 表达式只会被执行一次
  static std::once_flag flag;
  std::call_once(flag, [&interface, &refineMeshKeepStep]() {
    // reset interface and refineMeshKeepStep field
    interface = dimensionedScalar("zero", interface.dimensions(), 0.0);
    refineMeshKeepStep = dimensionedScalar("zero", refineMeshKeepStep.dimensions(), 0.0);
  });
  forAll (mesh_.C(), cellI) {
    // 当前 cellI 是否位于任意一个颗粒中
    bool cellInParticle = false;
    bool cellFirstEntryRefineMeshKeepStep = false;
    // 网格中心
    Foam::vector cellPos = mesh_.C()[cellI];
    for (int index = 0; index < numberOfParticles(); ++index) {
      if (checkPeriodicCells()) {
        FatalError << "Error: not support periodic check!" << abort(FatalError);
      }
      // 判断网格中心是否在 index 颗粒中
      double value = voidFractionM().pointInParticle(index, cellPos, refineMeshSkin());
      if (0 == refineMeshKeepInterval()) {
        if (value <= 0.0) {
          interface[cellI] = std::max(interface[cellI], value + 1);
        }
      } else { // valid refineMeshKeepInterval
        if (value <= 0.0) {
          // cellI 位于 index 颗粒内部
          cellInParticle = true;
          // 如果 cellI 网格位于任何一个颗粒内部，则重新设置 refineMeshKeepStep
          refineMeshKeepStep[cellI] = refineMeshKeepInterval();
          // 设置 interface
          interface[cellI] = std::max(interface[cellI], value + 1);
        } else {
          if (false == cellInParticle) {
            // cellI 目前不在任何颗粒中
            if (refineMeshKeepStep[cellI] > Foam::SMALL && false == cellFirstEntryRefineMeshKeepStep) {
              // refineMeshKeepStep[cellI] > 0.0，则保持 interFace 值
              cellFirstEntryRefineMeshKeepStep = true; // 确保对每一个 cellI 只执行一次
              refineMeshKeepStep[cellI] -= 1.0;
            } else {
              // 设置 interFace 为 0.0
              interface[cellI] = 0.0;
            }
          }
        }
      }
    } // end of all particles
  } // end of all cell in current processor
}

} // namespace Foam
