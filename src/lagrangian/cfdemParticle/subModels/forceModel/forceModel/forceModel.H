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
    forceModel

SourceFiles
    forceModel.C

\*---------------------------------------------------------------------------*/

#ifndef forceModel_H
#define forceModel_H

#include "fvCFD.H"
#include "cfdemCloud.H"
#include "probeModel.H"
#include "forceSubModel.H"
#include "interpolationCellPointFace.H"

namespace Foam {

class forceModel {

protected:

  const dictionary& dict_;

  cfdemCloud& particleCloud_;

  //! @brief 颗粒隐式力的总和 [N]
  mutable volVectorField impParticleForces_;

  //! @brief 颗粒显式力的总和 [N]
  mutable volVectorField expParticleForces_;

  //! @brief cfdemCloud 类型名称("A" or "B" or "Bfull" or "off")
  const word modelType_;

  //! @brief 是否激活探针
  bool probeIt_;

  //! @brief 是否需要方向矢量(orientation vector)
  bool requiresEx_;

  //! @brief 是否需要颗粒形状信息
  bool requiresShape_;

  bool requiresQuaternion_;

  bool requiresSuperquadric_;

  //! @brief 是否会改变转动信息(rotation information)
  bool pullPushRotation_;

  //! @brief 是否会改变动量交换场 Ksl
  bool implicitDrag_;

  //! @brief 是否会改变各向异性动量交换场 Ksl
  bool implicitAnisotropicDrag_;

  //! @brief 是否会改变转动动量交换场 Ksl
  bool implicitRotation_;

  //! @brief 保存当前 force model 的所有 forceSubModel 的名称
  // @note 每一类型的力都可以指定多个种类的 forceSubModel, 但是目前只有一种 forceSubModel 的实现, 即 ImEx
  wordList forceSubModels_;

  //! @brief 保存当前 force model 的所有 forceSubModel 的指针
  // @note 每一类型的力都可以指定多个种类的 forceSubModel, 但是目前只有一种 forceSubModel 的实现, 即 ImEx
  autoPtr<forceSubModel>* forceSubModel_;

  mutable autoPtr<interpolation<scalar> > GInterpolator_;

  mutable autoPtr<interpolation<scalar> > voidfractionInterpolator_;

  mutable autoPtr<interpolation<vector> > UInterpolator_;

  mutable autoPtr<interpolation<vector> > vorticityInterpolator_;

  mutable autoPtr<interpolation<vector> > gradPInterpolator_;

  mutable autoPtr<interpolation<vector> > gradUInterpolator_;

  mutable autoPtr<interpolation<vector> > gradVoidfractionInterpolator_;

  mutable autoPtr<interpolation<vector> > Up1Interpolator_;

  mutable autoPtr<interpolation<vector> > Up2Interpolator_;

  mutable autoPtr<interpolation<scalar> > dSauterInterpolator_;

  mutable autoPtr<interpolation<scalar> > phiP1Interpolator_;

  mutable autoPtr<interpolation<scalar> > phiP2Interpolator_;

  mutable autoPtr<interpolation<scalar> > alphaInterpolator_;

  mutable autoPtr<interpolation<vector> > gradAlphaInterpolator_;

  mutable autoPtr<interpolation<scalar> > TInterpolator_;

  mutable autoPtr<interpolation<vector> > UsInterpolator_;

  mutable autoPtr<interpolation<scalar> > fluidScalarFieldInterpolator_;

  mutable autoPtr<interpolation<vector> > gradPsolidInterpolator_;

  mutable autoPtr<interpolation<symmTensor> > shearRateInterpolator_;

  mutable autoPtr<interpolation<vector> > DDtUInterpolator_;

  mutable autoPtr<interpolation<vector> > divTauInterpolator_;

public:

  //! @brief Runtime type information
  TypeName("forceModel");

  //! @brief Declare runtime constructor selection table
  declareRunTimeSelectionTable(
    autoPtr,
    forceModel,
    dictionary,
    (
      const dictionary& dict,
      cfdemCloud& sm
    ),
    (dict,sm)
  );

  //! @brief Constructors
  forceModel(const dictionary& dict, cfdemCloud& sm);

  //! @brief Destructor
  virtual ~forceModel();

  //! @brief Selector
  static autoPtr<forceModel> New(const dictionary& dict,
                                 cfdemCloud& sm,
                                 word forceType);

  // Member Functions

  virtual void setForce() const = 0;

  virtual void setMixForce(const std::vector<double>& dimensionRatios) const {
    FatalError << "forceModel:setMixForce(): using base class function, please use derived class function\n" << abort(FatalError);
  }

  void applyDebugSettings(bool) const;

  //! @brief 对于不包含任何颗粒的网格使用显式力耦合
  void treatVoidCells() const;

  //! @brief 重新划分 Implicit / Explicit force
  void repartitionImExForces() const;

  //! @brief 设置 force sub model
  //! @brief 比如对于 DiFeliceDrag force model, setForceSubModels 会从 DiFeliceDragProps 中寻找并设置所有的 forceSubModels 字段
  //        如果 DiFeliceDragProps 中没有指定 forceSubModels 字段, 则会默认设置为 ImEx
  // @note 每一类型的力都可以指定多个种类的 forceSubModel, 但是目前只有一种 forceSubModel 的实现, 即 ImEx
  void setForceSubModels(dictionary& dict);

  //! @brief 添加由于多相作用而产生的湍流
  virtual void multiphaseTurbulence(volScalarField& field, bool) const { field *= 0.0; }

  //! @brief Function for simple explicit treatment of coupling terms, only for temperature
  virtual void manipulateScalarField(volScalarField& field) const { field *= 0.0; }

  //! @brief Function for explicit or implicit treatment of coupling terms, for heat and species balance equations
  virtual void manipulateScalarField(volScalarField& field, volScalarField& fieldImpl, int speciesID) const {
    field *= 0.0;
    fieldImpl *= 0.0;
  }

  //! @brief Function for return particle based data to DEM
  virtual void commToDEM() const {}

  virtual void setForcesOnParticle() const {}

  // Access Functions

  inline word modelType(){ return modelType_; }

  inline volVectorField& impParticleForces() const { return impParticleForces_; }

  inline volVectorField& expParticleForces() const { return expParticleForces_; }

  inline double** impForces() const { return particleCloud_.impForces_; }

  inline double** expForces() const { return particleCloud_.expForces_; }

  inline double** DEMForces() const { return particleCloud_.DEMForces_; }
  
  inline double** DEMTorques() const { return particleCloud_.DEMTorques(); }

  inline double** Cds() const { return particleCloud_.Cds_; }

  inline double** CdsExtra() const { return particleCloud_.CdsExtra(); }

  inline double** CdsRotation() const { return particleCloud_.CdsRotation(); }

  inline double** omegaFluid() const { return particleCloud_.omegaFluid(); }

  inline double** fluidVel() const { return particleCloud_.fluidVel_; }

  virtual inline bool& requiresEx() { return requiresEx_; }

  virtual inline bool& requiresQuaternion() { return requiresQuaternion_; }

  virtual inline bool& requiresSuperquadric() { return requiresSuperquadric_; }

  virtual inline bool& requiresShape() { return requiresShape_; }

  virtual inline bool& pullPushRotation() {return pullPushRotation_; }

  virtual inline bool& implicitDrag() {return implicitDrag_; }
  
  virtual inline bool& implicitAnisotropicDrag() {return implicitAnisotropicDrag_; }
  
  virtual inline bool& implicitRotation() {return implicitRotation_; }

  inline const wordList& forceSubModels() { return forceSubModels_; }

  inline const forceSubModel& forceSubM(int i) const { return forceSubModel_[i]; }

  inline int nrForceSubModels() const { return forceSubModels_.size(); }

};  // End of class forceModel

}  // End of namespace Foam

#endif  // forceModel_H
