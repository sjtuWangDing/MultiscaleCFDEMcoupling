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

    cloud class managing DEM data for CFD-DEM coupling

Class
    Foam::cfdemCloud

\*---------------------------------------------------------------------------*/

#ifndef __CFDEM_CLOUD_H__
#define __CFDEM_CLOUD_H__

#include <vector>
#include "OFversion.H" // choose version
#include "fvCFD.H"
#include "IFstream.H"

#if defined(version21) || defined(version16ext)
  #include "turbulenceModel.H"
#elif defined(version15)
  #include "RASModel.H"
#endif

#ifndef __MIXCLOUD__
#define __MIXCLOUD__ 1
#endif

namespace Foam {

// forward declarations
class forceModel;
class locateModel;
class voidFractionModel;
class dataExchangeModel;
class IOModel;
class probeModel;
class registryModel;
class averagingModel;
class clockModel;
class smoothingModel;
class momCoupleModel;
class meshMotionModel;
class liggghtsCommandModel;
class CouplingProperties;
class ParticleCloud;

class cfdemCloud {

public:

  //! @brief Construct from mesh and a list of particles
  cfdemCloud(const fvMesh& mesh);

  //! @brief Destructor
  virtual ~cfdemCloud();

protected:
  const fvMesh& mesh_;

  //! @note 在当前类中一定要最先声明 couplingPropertiesDict_ 和 liggghtsCommandsDict_
  IOdictionary couplingPropertiesDict_;

  IOdictionary liggghtsCommandsDict_;

  autoPtr<forceModel>* forceModel_;

  autoPtr<locateModel> locateModel_;

  autoPtr<momCoupleModel> momCoupleModel_;

  autoPtr<dataExchangeModel> dataExchangeModel_;

  autoPtr<IOModel> IOModel_;

  autoPtr<probeModel> probeModel_;

  autoPtr<registryModel> registryModel_;

  autoPtr<voidFractionModel> voidFractionModel_;

  autoPtr<averagingModel> averagingModel_;

  autoPtr<clockModel> clockModel_;

  autoPtr<smoothingModel> smoothingModel_;

  autoPtr<meshMotionModel> meshMotionModel_;

  autoPtr<liggghtsCommandModel> liggghtsCommand_;

  autoPtr<CouplingProperties> cProps_;

  autoPtr<ParticleCloud> pCloud_;

  bool writeTimePassed_;

  bool arraysReallocated_;

  /*!
   * @brief 是否隐式计算颗粒所受到的阻力
   *        true - 在每个耦合时间步，流体的速度和阻力系数都被传递到 DEM 中，从而在每个 DEM 时间步中，
   *               使用上一个耦合时间步中的阻力系数和流体速度，与当前颗粒速度一起计算颗粒受到的阻力
   *        false - 在每个耦合时间步中，流体对颗粒的阻力被传递到 DEM 中，并且在接下来的 DEM 时间步中，
   *                这个力保持不变，直到下个耦合时间步
   * @note default false
   */
  bool impDEMdrag_;

  /*!
   * @brief 在每个 DEM 时间步中，都将颗粒受到的力累计起来，然后在耦合时间步中传递给 CFD 计算，
   *        只有当 impDEMdrag_ = true 时才为真
   * @note default false
   */
  bool impDEMdragAcc_;

  //! @brief 显式、隐式分裂系数
  scalar impExpSplitFactor_;

  // volScalarField ddtVoidfraction_;

  //! @brief Variable used to de-activate mirroring across periodic boundary conditions.
  Switch checkPeriodicCells_;

  /*!
   * @brief De-activation and tolerance variables, if set to (for a specific direction),
   *        the periodic check will NOT be done. Important for probing ambient points.
   *        Only read-in in case checkPeriodicCells is active.
   * @note default = (1, 1, 1), i.e., periodic checks will be done
   */
  vector wallPeriodicityCheckRange_;

  scalar wallPeriodicityCheckTolerance_;

  bool meshHasUpdated_;

// #if defined(version24Dev)
//   const turbulenceModel& turbulence_;
// #elif defined(version21) || defined(version16ext)
//   #ifdef compre
//     const compressible::turbulenceModel& turbulence_;
//   #else
//     const incompressible::turbulenceModel& turbulence_;
//   #endif
// #elif defined(version15)
//   const incompressible::RASModel& turbulence_;
// #endif

  //! @brief Multiphase Turbulence (e.g., slip-induced turbulence)
  // volScalarField turbulenceMultiphase_;
};

} // namespace Foam

#endif // __CFDEM_CLOUD_H__
