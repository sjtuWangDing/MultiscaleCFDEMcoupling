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
  cfdemCloud class managing DEM data for CFD-DEM coupling

Class
  Foam::cfdemCloud
\*---------------------------------------------------------------------------*/

#ifndef __CFDEM_CLOUD_H__
#define __CFDEM_CLOUD_H__

#include <memory>
#include <vector>
#include <unordered_map>
#include <unistd.h> // usleep
#include "mpi.h" // MPI

#include "fvCFD.H"
#include "IFstream.H"
#include "error.h"

#include "base/tensor/tensor.H"
#include "cloud/OFVersion.H" // choose version
#include "cloud/couplingProperties.H"
#include "cloud/particleCloud.H"

#if defined(version21) || defined(version16ext)
  #include "turbulenceModel.H"
#elif defined(version15)
  #include "RASModel.H"
#endif

#ifndef USLEEP
#define USLEEP(times) \
  usleep((int)(times * 1000000)); \
  MPI_Barrier(MPI_COMM_WORLD);
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

class cfdemCloud {

public:

  //! @brief Construct from mesh and a list of particles
  cfdemCloud(const fvMesh& mesh);

  //! @brief Destructor
  virtual ~cfdemCloud();

  /*!
   * \brief 更新函数
   * \note used for cfdemSolverPiso
   * \param VoidF  <[in, out] 小颗粒空隙率场
   * \param Us     <[in, out] 局部平均小颗粒速度场
   * \param U      <[in] 流体速度场
   */
  virtual bool evolve(volScalarField& VoidF,
                      volVectorField& Us,
                      volVectorField& U);

  // interface of cfdemCloud
  inline const fvMesh& mesh() const { return mesh_; }
  inline ParticleCloud& pCloud() { return pCloud_; }
  inline const CouplingProperties& cProps() const { return cProps_; }
  inline const IOdictionary& couplingPropertiesDict() const { return couplingPropertiesDict_; }
  inline const IOdictionary& liggghtsCommandsDict() const { return liggghtsCommandsDict_; }

  const std::vector<std::shared_ptr<liggghtsCommandModel>>& liggghtsCommandModels() const;
  const std::vector<std::shared_ptr<forceModel>>& forceModels() const;
  const std::vector<std::shared_ptr<momCoupleModel>>& momCoupleModels() const;

  std::shared_ptr<forceModel> forceM(int index) const;
  const dataExchangeModel& dataExchangeM() const;
  const averagingModel& averagingM() const;
  const voidFractionModel& voidFractionM() const;
  dataExchangeModel& dataExchangeM();
  averagingModel& averagingM();
  voidFractionModel& voidFractionM();

  inline bool writeTimePassed() const { return writeTimePassed_; }

  // interface of CouplingProperties
  inline bool verbose() const { return cProps_.verbose(); }
  inline const word& modelType() const { return cProps_.modelType(); }
  inline const wordList& forceModelList() const { return cProps_.forceModelList(); }
  inline const wordList& momCoupleModelList() const { return cProps_.momCoupleModelList(); }
  inline const wordList& liggghtsCommandModelList() const { return cProps_.liggghtsCommandModelList(); }
  inline const word& turbulenceModelType() const { return cProps_.turbulenceModelType(); }
  inline const word& useDDTvoidfraction() const { return cProps_.useDDTvoidfraction(); }
  inline scalar impExpSplitFactor() const { return cProps_.impExpSplitFactor(); }
  inline bool solveFlow() const { return cProps_.solveFlow(); }
  inline bool treatVoidCellsAsExplicitForce() const { return cProps_.treatVoidCellsAsExplicitForce(); }
  inline bool debug() const { return cProps_.debug(); }
  inline bool ignore() const { return cProps_.ignore(); }
  inline bool allowAdjustTimeStep() const { return cProps_.allowAdjustTimeStep(); }
  inline bool allowCFDsubTimeStep() const { return cProps_.allowCFDsubTimeStep(); }
  inline label couplingInterval() const { return cProps_.couplingInterval(); }

  // interface of particleCloud
  inline int numberOfParticles() const { return pCloud_.numberOfParticles(); }

protected:

  const fvMesh& mesh_;

  //! @note 在当前类中一定要最先声明 couplingPropertiesDict_ 和 liggghtsCommandsDict_
  IOdictionary couplingPropertiesDict_;

  IOdictionary liggghtsCommandsDict_;

  CouplingProperties cProps_;

  ParticleCloud pCloud_;

  //! @note forceModel, liggghtsCommandModel, momCoupleModel using std::vector as container
  std::vector<std::shared_ptr<liggghtsCommandModel>> liggghtsCommandModels_;

  std::vector<std::shared_ptr<forceModel>> forceModels_;

  std::vector<std::shared_ptr<momCoupleModel>> momCoupleModels_;

  autoPtr<dataExchangeModel> dataExchangeModel_;

  autoPtr<averagingModel> averagingModel_;

  autoPtr<voidFractionModel> voidFractionModel_;

  autoPtr<locateModel> locateModel_;

  autoPtr<IOModel> IOModel_;

  autoPtr<probeModel> probeModel_;

  autoPtr<registryModel> registryModel_;

  autoPtr<clockModel> clockModel_;

  autoPtr<smoothingModel> smoothingModel_;

  autoPtr<meshMotionModel> meshMotionModel_;

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
