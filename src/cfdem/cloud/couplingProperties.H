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
    Foam::CouplingProperties
\*---------------------------------------------------------------------------*/

#ifndef __COUPLING_PROPERTIES_H__
#define __COUPLING_PROPERTIES_H__

#include "cfdemCloud.H"

namespace Foam {

class CouplingProperties {

public:

  CouplingProperties(const fvMesh& mesh,
                     const IOdictionary& couplingPropertiesDict,
                     const IOdictionary& liggghtsCommandsDict);

  // interfaces - getter
  inline bool verbose() const { return verbose_; }
  inline const word& modelType() const { return modelType_; }
  inline const wordList& forceModels() const { return forceModels_; }
  inline const wordList& momCoupleModels() const { return momCoupleModels_; }
  inline const wordList& liggghtsCommandModels() const { return liggghtsCommandModels_; }
  inline const word& turbulenceModelType() const { return turbulenceModelType_; }

protected:
  const fvMesh& mesh_;

  //! @note 在当前类中一定要最先声明 couplingPropertiesDict_ 和 liggghtsCommandsDict_
  const IOdictionary& couplingPropertiesDict_;

  const IOdictionary& liggghtsCommandsDict_;

  //! @brief 是否打印多余调试信息
  bool verbose_;

  //! @brief 模型类型(A, B, Bfull, none)
  const word modelType_;

  //! @brief dict 中指定的所有 force model 的名称
  const wordList forceModels_;

  //! @brief dict 中指定的所有 mom couple model 的名称
  const wordList momCoupleModels_;

  //! @brief dict 中指定的所有 liggghts command model 的名称
  const wordList liggghtsCommandModels_;

  //! @brief dict 中指定的湍流模型的名称
  const word turbulenceModelType_;

  //! @brief 显式、隐式分裂系数
  scalar impExpSplitFactor_;

  //! @brief 是否需要求解流体
  Switch solveFlow_;

  //! @brief Treat the coupling force in cells with no Us data explicitly.
  bool treatVoidCellsAsExplicitForce_;

  bool debug_;

  bool ignore_;

  bool allowAdjustTimeStep_;

  word useDDtVoidfraction_;

#if __MIXCLOUD__

public:

  //! @brief 颗粒尺度系数
  static double fineParticleRatio_;

  static double coarseParticleRatio_;

protected:

  //! @brief 是否用于求解器 cfdemSolverIB
  bool usedForSolverIB_;

  //! @brief 是否用于求解器 cfdemSolverPiso
  bool usedForSolverPiso_;

  //! @brief 是否使用 dynamic refine mesh
  bool useDynamicRefineMesh_;

  //! @brief true if particle is fixed
  bool fixedParticle_;

  //! @brief 来流速度
  vector flowVelocity_;

public:

  inline bool checkCoarseParticle(const double& ratio) const {
    return ratio > 0 && ratio < coarseParticleRatio_;
  }
  inline bool checkMiddleParticle(const double& ratio) const {
    return ratio > 0 && (coarseParticleRatio_ <= ratio && ratio < fineParticleRatio_);
  }
  inline bool checkFineParticle(const double& ratio) const {
    return ratio > 0 && ratio >= fineParticleRatio_;
  }
  inline bool checkFineAndMiddleParticle(const double& ratio) const {
    return !checkCoarseParticle(ratio);
  }
  /*!
   * \brief 是否需设置 field for coarse particle
   * \param index            <[in] 颗粒索引
   * \param dimensionRatios  <[in] 颗粒尺度比值
   */
  inline bool needSetFieldForCoarseParticle(const int& index,
                                            const std::vector<double> dimensionRatios) const {
    if (usedForSolverIB_ || useDynamicRefineMesh_) {
      return true;
    } else if (!dimensionRatios.empty() && !usedForSolverPiso_) {
      return index < static_cast<int>(dimensionRatios.size()) ? checkCoarseParticle(dimensionRatios[index]) : false;
    } else {
      return false;
    }
  }
  inline bool usedForSolverIB() const { return usedForSolverIB_; }
  inline bool usedForSolverPiso() const { return usedForSolverPiso_; }
  inline bool useDynamicRefineMesh() const { return useDynamicRefineMesh_; }
  inline bool fixedParticle() const { return fixedParticle_; }
  inline const vector& flowVelocity() const { return flowVelocity_; }

#endif // __MIXCLOUD__

};

} // namespace Foam

#endif // __COUPLING_PROPERTIES_H__
