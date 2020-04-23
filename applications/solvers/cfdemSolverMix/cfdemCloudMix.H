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

    cloud class managing DEM data for CFD-DEM coupling and IB representation

Class
    Foam::cfdemCloudMix derived from cfdemCloud

SourceFiles
    cfdemCloudMix.C

\*---------------------------------------------------------------------------*/

#ifndef cfdemCloudMix_H
#define cfdemCloudMix_H

#include "cfdemCloud.H"
#include <unistd.h>
#include <functional>

namespace Foam {

class cfdemCloudMix: public cfdemCloud {

protected:

  // @brief 颗粒角速度
  mutable double** angularVelocities_;

  // @brief 流体对颗粒的力矩
  mutable double** DEMTorques_;

  // @brief 是否跳过拉格朗日映射到欧拉场
  mutable bool skipLagrangeToEulerMapping_;

  // @brief 是否在指定时间步后 Skip Lagrange To Euler Mapping.
  mutable bool skipAfter_;

  // @brief 指定在多少时间步后 Skip Lagrange To Euler Mapping.
  mutable int  timeStepsToSkip_;

  mutable bool calculateTortuosity_;

  // @brief Mean flow direction
  vector flowDir_;

  volScalarField* frontMeshRefineField_;

  bool frontMeshRefine_;

  label pRefCell_;

  scalar pRefValue_;

  mutable bool haveEvolvedOnce_;

  // @brief 颗粒体积分数
  mutable double** volumefractions_;

  // @brief 颗粒网格尺寸比
  std::vector<double> dimensionRatios_;

  // @brief 全局颗粒网格尺寸比
  std::vector<double> globalDimensionRatios_;

  // @brief 颗粒覆盖当前求解器上网格的总个数
  std::vector<double> sumCellsNumbers_;

  // @brief 颗粒覆盖当前求解器上网格的总体积
  std::vector<double> sumCellsVolumes_;

  // // @brief fine 颗粒索引
  // std::vector<int> fineParticleIndexs_;

  // // @brief middle 颗粒索引
  // std::vector<int> middleParticleIndexs_;

  // // @brief coarse 颗粒索引
  // std::vector<int> coarseParticleIndexs_;

  // @brief 颗粒影响的网格集合
  // @note 用于计算 drag force for middle particle
  std::vector<std::vector<int> > affectedCellIDs;

private:

  // Functions for debug
  void doCoupleDebug();

  // Functions for debug
  void findCellDebug();

public:

  /*! \brief Constructors */
  cfdemCloudMix(const fvMesh& mesh);

  /*! \brief Destructor */
  ~cfdemCloudMix();

  virtual void giveDEMdata();

  virtual void getDEMdata();

  virtual bool reAllocArrays() const;

  /*! \brief force number of particles during reallocation */
  virtual bool reAllocArrays(int nP, bool forceRealloc) const;

  /*! \brief 强制重新分配内存 */
  void mixForceReAllocArrays();

  void mixResetFieldKernel();

  void mixCouplingKernel();

  /*!
   * \brief 更新函数
   * \note used for cfdemSolverPiso
   * \param alphaVoidfraction  <[in, out] 小颗粒空隙率场
   * \param Us                 <[in, out] 局部平均小颗粒速度场
   * \param U                  <[in] 流体速度场
   */
  virtual bool evolve(volScalarField& alphaVoidfraction,
                      volVectorField& Us,
                      volVectorField& U);

  /*!
   * \brief 更新函数
   * \note used for cfdemSolverIB
   * \param alpha      <[in, out] 大颗粒体积分数场
   * \param interFace  <[in, out] 界面场
   */
  bool evolve(volScalarField& alpha,
              volScalarField& interFace);

  /*!
   * \brief 更新函数
   * \note used for cfdemSolverMix
   * \param alphaVoidfraction    <[in, out] 小颗粒空隙率场
   * \param alphaVolumefraction  <[in, out] 大颗粒体积分数场
   * \param Us                   <[in, out] 局部平均小颗粒速度场
   * \param U                    <[in] 流体速度场
   * \param interFace            <[in, out] 界面场
   */
  bool evolve(volScalarField& alphaVoidfraction,
              volScalarField& alphaVolumefraction,
              volVectorField& Us,
              volVectorField& U,
              volScalarField& interFace);

  void calcVelocityCorrection(volScalarField& p,
                              volVectorField& U,
                              volScalarField& phiIB,
                              volScalarField& voidfraction);

  // @brief 设置 coarse particle 中的流体速度
  void setParticleVelocity(volVectorField& U);

  // @brief 确定颗粒周围细化网格的区域(每个方向的尺寸都是颗粒尺寸的两倍)
  void setInterFace(volScalarField& interFace);

  // Access Functions
  inline double** angularVelocities() const { return angularVelocities_; }

  inline double** DEMTorques() const { return DEMTorques_; }

  inline vector angularVelocity(int index) {
    vector vel;
    for(int i = 0; i < 3; i++) { vel[i] = angularVelocities_[index][i]; }
    return vel;
  }

};  // End of class cfdemCloudMix

}  // End of namespace Foam

#endif  // cfdemCloudMix_H
