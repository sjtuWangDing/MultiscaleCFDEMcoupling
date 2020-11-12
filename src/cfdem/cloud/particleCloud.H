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
#include "base/tenser/tenser.H"

namespace Foam {

//! @brief The basic class for particles
class ParticleCloud {

public:

  ParticleCloud(int numberOfParticle = -1);

  // interfaces - getter
  inline int numberOfParticles() const { return numberOfParticles_; }
  inline bool numberOfParticlesChanged() const { return numberOfParticlesChanged_; }
  inline double** positions() const { return positions_; }
  inline double** velocities() const { return velocities_; }
  inline double** fluidVelocites() const { return fluidVelocites_; }
  inline double** cds() const { return cds_; }
  inline double** fAcc() const { return fAcc_; }
  inline double** impForces() const { return impForces_; }
  inline double** expForces() const { return expForces_; }
  inline double** demForces() const { return demForces_; }
  inline double** radii() const { return radii_; }
  inline double** volumes() const { return volumes_; }
  inline double** voidfractions() const { return voidfractions_; }
  inline double** cellIds() const { return cellIds_; }
  inline double** particleWeights() const { return particleWeights_; }
  inline double** particleVolumes() const { return particleVolumes_; }

  // interface - setter
  inline void setNumberOfParticles(int number) {
    if (number <= 0) {
      FatalError << "number of particle is <= 0" << abort(FatalError);
    }
    numberOfParticles_ = number;
  }

protected:

  //! @brief 颗粒总数
  int numberOfParticles_;

  //! @brief 颗粒总数是否发生变化
  bool numberOfParticlesChanged_;

  //! @brief 颗粒位置矢量
  //! @note positions_[index][0 ~ 2]
  double** positions_;

  //! @brief 颗粒速度
  //! @note velocities_[index][0 ~ 2]
  double** velocities_;

  //! @brief 小颗粒的中心处的流体速度
  //         该速度用于传递给 DEM 求解器, 则 DEM 求解器使用颗粒中心处的流体速度与阻力系数一起计算颗粒受到的阻力
  //! @note fluidVelocites_[index][0 ~ 2]
  double** fluidVelocites_;

  //! @brief 颗粒阻力系数
  //         该阻力系数用于传递给 DEM 求解器, 则 DEM 求解器使用颗粒中心处的流体速度与阻力系数一起计算颗粒受到的阻力
  //! @note cds_[index][0]
  double** cds_;

  //! @brief 在一个耦合时间间隔中，小颗粒所受的阻力的累计，将在下个耦合步中，传递给流场
  //! @note fAcc_[index][0 ~ 2] ??????????
  double** fAcc_;

  //! @brief 颗粒对流体的隐式作用力
  //! @note impForces_[index][0 ~ 2]
  double** impForces_;

  //! @brief 颗粒对流体的显式作用力
  //! @note expForces_[index][0 ~ 2]
  double** expForces_;

  //! @brief 流体对颗粒的总作用力
  //! @note demForces_[index][0 ~ 2] ??????????
  double** demForces_;

  //! @brief 颗粒半径
  //! @note radii_[index][0]
  double** radii_;

  //! @brief 颗粒体积
  //! @note volumes_[index][0]
  double** volumes_;

  //! @brief 空隙率
  //! @note voidfractions_[index][subcell]
  double** voidfractions_;

  //! @brief 颗粒覆盖的所有网格的编号
  //! @note cellIds_[index][subcell]
  double** cellIds_;

  //! @brief 颗粒对所覆盖网格的影响系数
  //         如果使用 divided 空隙率模型，则对颗粒覆盖的某个网格 subcell，
  //         只要有一个颗粒标志点在网格中，particleWeights_[index][subcell] += 1.0 / 29.0
  //! @note particleWeights_[index][subcell]
  double** particleWeights_;

  //! @brief 颗粒对所覆盖网格的覆盖体积
  //         如果使用 divided 空隙率模型，则对颗粒覆盖的某个网格 subcell，
  //         只要有一个颗粒标志点在网格中，particleVolumes_[index][subcell] += (1.0 / 29.0) * 颗粒体积
  //! @note particleVolumes_[index][subcell]
  double** particleVolumes_;
};

} // namespace Foam

#endif // __PARTICLE_CLOUD_H__
