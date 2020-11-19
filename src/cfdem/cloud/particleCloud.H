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
  Foam::ParticleCloud
\*---------------------------------------------------------------------------*/

#ifndef __PARTICLE_CLOUD_H__
#define __PARTICLE_CLOUD_H__

#include "error.H"
#include "base/tensor/tensor.H"
#include "cloud/cfdemBase.H"

namespace Foam {

//! \brief The basic class for particles
class ParticleCloud {

public:

  ParticleCloud(int numberOfParticle = -1);

  // interfaces - getter
  inline int numberOfParticles() const { return numberOfParticles_; }
  inline bool numberOfParticlesChanged() const { return numberOfParticlesChanged_; }
  inline const base::CITensor1& particleOverMeshNumber() const { return particleOverMeshNumber_; }
  inline const base::CDTensor1& radii() const { return radii_; }
  inline const base::CDTensor1& volumes() const { return volumes_; }
  inline const base::CDTensor1& cds() const { return cds_; }
  inline const base::CITensor2& cellIDs() const { return cellIDs_; }
  inline const base::CDTensor2& positions() const { return positions_; }
  inline const base::CDTensor2& velocities() const { return velocities_; }
  inline const base::CDTensor2& fluidVel() const { return fluidVel_; }
  inline const base::CDTensor2& fAcc() const { return fAcc_; }
  inline const base::CDTensor2& impForces() const { return impForces_; }
  inline const base::CDTensor2& expForces() const { return expForces_; }
  inline const base::CDTensor2& DEMForces() const { return DEMForces_; }
  inline const base::CDTensor2& DEMTorques() const { return DEMTorques_; }
  inline const base::CDTensor2& voidfractions() const { return voidfractions_; }
  inline const base::CDTensor2& particleWeights() const { return particleWeights_; }
  inline const base::CDTensor2& particleVolumes() const { return particleVolumes_; }

  // interface - setter
  inline void setNumberOfParticles(int number) {
    if (number <= 0) {
      FatalError << "number of particle is <= 0" << abort(FatalError);
    }
    numberOfParticles_ = number;
  }

private:

  //! \brief 颗粒总数
  int numberOfParticles_;

  //! \brief 颗粒总数是否发生变化
  bool numberOfParticlesChanged_;

#if CFDEM_USE_TENSOR
  // cpu int tensor 1D
  base::CITensor1 particleOverMeshNumber_;
  // cpu double tensor 1D
  base::CDTensor1 radii_;
  base::CDTensor1 volumes_;
  base::CDTensor1 cds_;
  // cpu int tensor 2D
  base::CITensor2 cellIDs_;
  // cpu double tensor 2D
  base::CDTensor2 positions_;
  base::CDTensor2 velocities_;
  base::CDTensor2 fluidVel_;
  base::CDTensor2 fAcc_;
  base::CDTensor2 impForces_;
  base::CDTensor2 expForces_;
  base::CDTensor2 DEMForces_;
  base::CDTensor2 DEMTorques_;
  base::CDTensor2 voidfractions_;
  base::CDTensor2 particleWeights_;
  base::CDTensor2 particleVolumes_;
#else
  //! \brief 颗粒半径
  //! \note radii_[index][0]
  double** radii_;

  //! \brief 颗粒体积
  //! \note volumes_[index][0]
  double** volumes_;

  //! \brief 颗粒阻力系数
  //         当 forceSubModel implForceDEM 为 true 时候，该 cds_ 和 fluidVel_ 传递给 DEM 求解器， 
  //         DEM 求解器使用颗粒中心处的流体速度与阻力系数一起计算颗粒受到的阻力
  //! \note cds_[index][0]
  double** cds_;

  //! \brief 颗粒位置矢量
  //! \note positions_[index][0 ~ 2]
  double** positions_;

  //! \brief 颗粒速度
  //! \note velocities_[index][0 ~ 2]
  double** velocities_;

  //! \brief 小颗粒的中心处的流体速度
  //         当 forceSubModel implForceDEM 为 true 时候，该 cds_ 和 fluidVel_ 传递给 DEM 求解器， 
  //         DEM 求解器使用颗粒中心处的流体速度与阻力系数一起计算颗粒受到的阻力
  //! \note fluidVel_[index][0 ~ 2]
  double** fluidVel_;

  //! \brief 在一个耦合时间间隔中，小颗粒所受的阻力的累计，将在下个耦合步中，传递给流场
  //! \note fAcc_[index][0 ~ 2] ??????????
  double** fAcc_;

  //! \brief 颗粒对流体的隐式作用力
  //! \note impForces_[index][0 ~ 2]
  double** impForces_;

  //! \brief 颗粒对流体的显式作用力
  //! \note expForces_[index][0 ~ 2]
  double** expForces_;

  //! \brief 流体对颗粒的总作用力，如果不设置 cds_，则
  //! \note DEMForces_[index][0 ~ 2]
  double** DEMForces_;

  //! \brief 空隙率
  //! \note voidfractions_[index][subcell]
  double** voidfractions_;

  //! \brief 颗粒覆盖的所有网格的编号
  //! \note cellIDs_[index][subcell]
  double** cellIDs_;

  //! \brief 颗粒对所覆盖网格的影响系数
  //         如果使用 divided 空隙率模型，则对颗粒覆盖的某个网格 subcell，
  //         只要有一个颗粒标志点在网格中，particleWeights_[index][subcell] += 1.0 / 29.0
  //! \note particleWeights_[index][subcell]
  double** particleWeights_;

  //! \brief 颗粒对所覆盖网格的覆盖体积
  //         如果使用 divided 空隙率模型，则对颗粒覆盖的某个网格 subcell，
  //         只要有一个颗粒标志点在网格中，particleVolumes_[index][subcell] += (1.0 / 29.0) * 颗粒体积
  //! \note particleVolumes_[index][subcell]
  double** particleVolumes_;
#endif // CFDEM_CLOUD_USE_TENSOR
};

} // namespace Foam

#endif // __PARTICLE_CLOUD_H__
