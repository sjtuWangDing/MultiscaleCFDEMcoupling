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
    forceSubModel

SourceFiles
    forceSubModel.C

\*---------------------------------------------------------------------------*/

#ifndef forceSubModel_H
#define forceSubModel_H

#include "fvCFD.H"
#include "cfdemCloud.H"
#include "probeModel.H"

namespace Foam {

/*
forceSubModel 中指定的 11 个 switch:

(0) treatForceExplicit
  true - 在CFD的动量方程中，耦合力都为显式力
  false (default) - 在CFD的动量方程中，耦合力都为半隐式力
  
(1) treatForceDEM
  true - 仅在 DEM 中考虑耦合力，即仅考虑流体对颗粒的作用，不考虑颗粒对流体的反作用
  false (default) - 即 CFD 和 DEM 都考虑耦合力

(2) implForceDEM
  true - 在每个耦合时间步，流体的速度和阻力系数都被传递到 DEM 中，从而在每个 DEM 时间步中，使用上一个耦合时间步中的阻力系数和流体速度，与当前颗粒速度一起计算颗粒受到的阻力
  false (default) - 在每个耦合时间步中，流体对颗粒的阻力被传递到 DEM 中，并且在接下来的 DEM 时间步中，这个力保持不变，直到下个耦合时间步

(3) verbose
  true - 调试信息输出到屏幕
  false (default) - 调试信息不输出到屏幕

(4) interpolation
  true - 将欧拉场使用插值模型插值到拉格朗日场
  false (default) - 只使用网格中心的值

(5) useFilteredDragModel
  true or false (default) - 是否使用 Beetstra 阻力模型的粗网格版本(coarse grid version)，考虑网格尺度对颗粒所受阻力的影响

(6) useParcelSizeDependentFilteredDrag
  true or false (default) - 是否使用 Beetstra 阻力模型的粗网格版本(coarse grid version)，考虑颗粒尺度对颗粒所受阻力的影响，同时，强制使 useFilteredDragModel 为 true

(7) implForceDEMaccumulated
  note - 只能与 implForceDEM 一起使用
  true or false (default) - 是否在每个 DEM 时间步中将颗粒受到的力累计起来，然后在耦合时间步中传递给 CFD 计算

(8) scalarViscosity
  true or false (default) - 是否使用用户自定义动力粘度 nu 进行阻力计算，而 CFD 计算仍然使用 transport dict 指定的动力粘度

(9) verboseToDisk
  true or false (default) - 是否将调试信息输出到文件

(10) useCorrectedVoidage
  true or false (default) - 是否使用修正后的空隙率
  */

class forceSubModel {

protected:

  const dictionary& dict_;

  cfdemCloud& particleCloud_;

  forceModel& forceModel_;

  // \brief 在基类中被定义的 switch 个数，共有 11 个
  label nrDefaultSwitches_;

  // \brief switch 变量的名称链表，共有 11 个元素
  wordList switchesNameList_;

  // \brief 设置每个 switch 是否将会被搜索，共有 11 个 Switch 类型的元素
  // @usage 如果 switchesList_[i] 被设置为了 true，则就会在 dict 中搜索名为 switchesNameList_[i] 的 switch，如果在 dict 中，指定为 true，则将 switches_[i] 设置为 true，如果指定为 false，则将 switches_[i] 设置为 false
  mutable List<Switch> switchesList_;

  // \brief 记录每个 switch 的设置值，用于计算，共有 11 个 Switch 类型的元素
  mutable List<Switch> switches_;

  // \brief 用户自定义动力粘度，该粘度仅用于计算颗粒所受阻力
  mutable volScalarField nu_;

  // \brief 计算 div(Tau_)
  mutable volVectorField divTau_;

  // \brief 单位体积上的 IB 阻力
  // @usage IBDragPerV_[cellI]，cellI 为流体网格编号
  mutable volVectorField IBDragPerV_;

  word densityFieldName_;

  const volScalarField& rho_;

  // \brief 调试信息输出到文件的间隔，初值是 1
  mutable int verboseDiskIntervall_;

  // \brief 调试信息输出到文件的计数，初值为 0
  mutable int verboseDiskCounter_;

  scalar scaleDia_;

  scalar scaleDrag_;

  scalar scaleDH_;

public:

  // \brief Runtime type information
  TypeName("forceSubModel");

  // \brief Declare runtime constructor selection table
  declareRunTimeSelectionTable(
    autoPtr,
    forceSubModel,
    dictionary,
    (
      const dictionary& dict,
      cfdemCloud& sm,
      forceModel& fm
    ),
    (dict, sm, fm)
  );

  // \brief Constructors
  forceSubModel(const dictionary& dict, cfdemCloud& sm, forceModel& fm);

  // \brief Destructor
  virtual ~forceSubModel();

  // \brief Selector
  static autoPtr<forceSubModel> New(const dictionary& dict,
                                    cfdemCloud& sm,
                                    forceModel& fm,
                                    word forceType);

  // Member Functions

  // \brief 计算索引为 index 的颗粒的 scale 直径
  // \param d     <[in, out] 颗粒直径
  // \param index <[in] 颗粒索引
  virtual void scaleDia(scalar& d, int index = 0) const;

  virtual void scaleForce(vector& force, scalar& d, int index = 0) const;

  virtual void scaleCoeff(scalar& coeff, scalar& d, int index = 0) const;

  // \param index                  <[in] 颗粒索引
  // \param dragTot                <[in] 索引为 index 的颗粒受到的总阻力
  // \param dragEx                 <[in] 索引为 index 的颗粒受到的显式阻力
  // \param Ufluid = vector::zero  <[in] 索引为 index 的颗粒中心处流体速度(可以指定是否使用插值模型计算)
  // \param scalar Cd = 0          <[in] 颗粒阻力系数
  void partToArray(const label& index,
                   const vector& dragTot,
                   const vector& dragEx,
                   const vector& Ufluid = vector::zero,
                   scalar Cd = scalar(0)) const;

  void readSwitches() const;

  void setSwitchesList(label i, bool v) const { switchesList_[i] = v; }

  void setSwitches(label i, Switch v) const { switches_[i] = v; }

  const volScalarField& nuField() const;

  const volScalarField& muField() const;

  const volScalarField& rhoField() const;

  const volVectorField& divTauField(const volVectorField& U) const;

  const volVectorField& IBDragPerV(const volVectorField& U, const volScalarField& p) const;

  // \brief function to set anisotropic drag components for DEM side only.
  void partToArrayAnisotropic(const label& index,
                              const vector& CdExtra,
                              const vector& dragEx = vector::zero) const {
    FatalError << "forceSubModel::partToArrayAnisotropic(): not implemented, please refer to forceSubModel.C in https://github.com/CFDEMproject/CFDEMcoupling-PUBLIC" << endl;
  }

  // \brief function to set anisotropic torque components for DEM side only.
  void partToArrayAnisotropicTorque(const label& index,
                                    const vector& CdTorque,
                                    const vector& torqueTotal = vector::zero) const {
    FatalError << "forceSubModel::partToArrayAnisotropicTorque(): not implemented, please refer to forceSubModel.C in https://github.com/CFDEMproject/CFDEMcoupling-PUBLIC" << endl;
  }

  virtual void constructorCalls(word typeName) const {}

  virtual void preParticleLoop(bool verbose = false) const {}

  virtual void postParticleLoop(bool verbose = false) const {}

  virtual void update(label particleI,
                      label cellI,
                      scalar& d,
                      scalar& scalToUpdate1,
                      scalar& scalToUpdate2,
                      bool verbose) const {}

  virtual void update(label particleI,
                      label cellI,
                      scalar& d,
                      vector& vecToUpdate1,
                      vector& vecToUpdate2,
                      scalar& scalToUpdate1,
                      scalar& scalToUpdate2,
                      bool verbose) const {}

  virtual void explicitCorr(vector& dragImplicit,
                            vector& dragExplicit,
                            scalar& dragCoefficient,
                            vector& Ufluid,
                            const vector& Ucell,
                            vector& Us,
                            const vector& UsCell,
                            bool verbose,
                            label index = 100) const {
    dragExplicit = vector::zero;
  }

  virtual void explicitCorrScalar(scalar& sourceKImplicit,
                                  scalar& sourceExplicit,
                                  scalar& areaTimesTransferCoefficient,
                                  scalar& fluidProperty,
                                  const scalar& fluidPropertyCell,
                                  scalar& particleProperty,
                                  bool verbose,
                                  label index = 100) const {
    sourceExplicit = areaTimesTransferCoefficient * (fluidProperty - particleProperty);
    sourceKImplicit = 0.0;
  }

  virtual void explicitLimit(vector& dragImplicit,
                             vector& dragExplicit,
                             scalar& d) const {
    dragExplicit = vector::zero;
  }

  virtual void verboseToDiskWrite(Field<scalar>& writeValues) const {}

  virtual void updateField(volScalarField* scalField, volVectorField* vecField) const {};

  virtual void getField(vector position, label cellI, scalar& scalVal, vector& vecVal) const {};

  virtual scalar calculateCorrectedVoidage(scalar voidfraction, scalar cellDpRatio) const { return voidfraction; }

  // Access Functions

  virtual word myType() const = 0;

  inline forceModel& myForceM() const { return forceModel_; }

  const List<Switch>& switches() const { return switches_; }

  const wordList& switchesNameList() const { return switchesNameList_; }

  inline bool treatForceExplicit() const { return switches_[0]; }

  inline bool treatForceDEM() const { return switches_[1]; }

  inline bool implForceDEM() const { return switches_[2]; }

  inline bool verbose() const { return switches_[3]; }

  inline bool interpolation() const { return switches_[4]; }

  inline bool useFilteredDragModel() const { return switches_[5]; }

  inline bool useParcelSizeDependentFilteredDrag() const { return switches_[6]; }

  inline bool implForceDEMaccumulated() const { return switches_[7]; }

  inline bool scalarViscosity() const { return switches_[8]; }

  inline bool verboseToDisk() const {
    verboseDiskCounter_ += 1;
    if (verboseDiskCounter_ >= verboseDiskIntervall_) {
      verboseDiskCounter_ = 0;
      return switches_[9];
    } else {
      return false;
    }
  }

  inline bool useCorrectedVoidage() const { return switches_[10]; }

  const label& nrDefaultSwitches() const { return nrDefaultSwitches_; }

  const scalar& scaleDia() const { return scaleDia_; }

  const scalar& scaleDrag() const { return scaleDia_; }

};  // End of class forceSubModel

}  // End of namespace Foam

#endif  // forceSubModel_H
