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
    averagingModel

SourceFiles
    averagingModel.C
\*---------------------------------------------------------------------------*/

#ifndef averagingModel_H
#define averagingModel_H

#include "fvCFD.H"
#include "cfdemCloud.H"
#include "dataExchangeModel.H"
#include "voidFractionModel.H"

namespace Foam {

//! @brief 局部平均模型
class averagingModel {

protected:

  const dictionary& dict_;

  cfdemCloud& particleCloud_;

  //! @brief 覆盖 cellI 网格的所有颗粒对其影响系数的总和
  //! @note UsWeightField_[cellI]
  mutable volScalarField UsWeightField_;

  //! @brief 上个时间步中局部平均颗粒速度场
  mutable volVectorField UsPrev_;

  //! @brief 下个时间步中局部平均颗粒速度场
  mutable volVectorField UsNext_;

public:

  //! @brief Runtime type information
  TypeName("averagingModel");

  //! @brief Declare runtime constructor selection table
  declareRunTimeSelectionTable(
    autoPtr,
    averagingModel,
    dictionary,
    (
      const dictionary& dict,
      cfdemCloud& sm
    ),
    (dict, sm)
  );

  //! @brief Constructors
  averagingModel(const dictionary& dict, cfdemCloud& sm);

  //! @brief Destructor
  virtual ~averagingModel();

  //! @brief Selector
  static autoPtr<averagingModel> New(const dictionary& dict, cfdemCloud& sm);

  void applyDebugSettings(bool debug) const;

  /*!
   * \brief 设置局部平均标量场(纯虚函数，由子类提供覆盖版本)
   * \param field                      <[in, out] 需要被局部平均化场
   * \param value                      <[in] 用于局部平均化的颗粒(lagrange)变量
   * \param weight                     <[in] 用于局部平均化的权重系数(lagrange)
   * \param weightField                <[in, out] 权重系数平均化场
   * \param mask
   * \param weight2                    <[in] 指定第二权重系数
   * \param weightWithWeight2 = false  <[in] 是否使用第二权重系数
   */
  virtual void setScalarAverage(volScalarField& field,
                                double**& value,
                                double**& weight,
                                volScalarField& weightField,
                                double** const& mask,
                                double** const& weight2,
                                bool weightWithWeight2 = false) const = 0;

  /*!
   * \brief 设置局部平均矢量场(纯虚函数，由子类提供覆盖版本)
   * \param field                      <[in, out] 需要被局部平均化场
   * \param value                      <[in] 用于局部平均化的颗粒(lagrange)变量
   * \param weight                     <[in] 用于局部平均化的权重系数(lagrange)
   * \param weightField                <[in, out] 权重系数平均化场
   * \param mask
   * \param weight2                    <[in] 指定第二权重系数
   * \param weightWithWeight2 = false  <[in] 是否使用第二权重系数
   */
  virtual void setVectorAverage(volVectorField& field,
                                double**& value,
                                double**& weight,
                                volScalarField& weightField,
                                double** const& mask,
                                double** const& weight2,
                                bool weightWithWeight2 = false) const = 0;

  /*!
   * \brief 设置局部平均矢量场
   * \param value              <[in] 用于局部平均化的颗粒(lagrange)变量
   * \param weight             <[in] 用于局部平均化的权重系数(lagrange)
   * \param dimensionRatios    <[in] 颗粒网格尺寸比
   * \param weightField        <[in, out] 权重系数平均化场
   * \param field              <[in, out] 需要被局部平均化场
   * \param mask
   */
  virtual void setMixVectorAverage(double **& value,
                                   double **& weight,
                                   std::vector<double>& dimensionRatios,
                                   volVectorField& field,
                                   volScalarField& weightField,
                                   double **const& mask) const {
    FatalError << "averagingModel::setMixVectorAverage() not implemented, please use dense model."
      << abort(FatalError);
  }

  void setVectorSum(volVectorField& field,
                    double**& value,
                    double**& weight,
                    double** const& mask) const;

  void setVectorSum(volVectorField& field,
                    double**& value,
                    double**& weight,
                    const std::vector<double>& dimensionRatios) const;

  void setScalarSum(volScalarField& field,
                    double**& value,
                    double**const& weight,
                    double**const& mask) const;

  // * * * * * * * * * * * * * *  Not Implemented * * * * * * * * * * * * * * * * //

  // @note If need, refer to averagingModel.C in https://github.com/CFDEMproject/CFDEMcoupling-PUBLIC
  void setVectorSumSimple(volVectorField& field,
                          double**& value,
                          double**& weight,
                          int nP) const {}

  // @note If need, refer to averagingModel.C in https://github.com/CFDEMproject/CFDEMcoupling-PUBLIC
  void undoVectorAverage(volVectorField& fieldPrev,
                         volVectorField& fieldNext,
                         volScalarField& weightField,
                         double** const& value,
                         double** const& weight,
                         double**const& mask,
                         bool single = false) const {}

  // @note If need, refer to averagingModel.C in https://github.com/CFDEMproject/CFDEMcoupling-PUBLIC
  void undoVectorSum(volVectorField& field,
                     double** const& value,
                     double** const& weight,
                     double**const& mask) const {}

  // @note If need, refer to averagingModel.C in https://github.com/CFDEMproject/CFDEMcoupling-PUBLIC
  void setDSauter(volScalarField& dSauter,
                  double**& weight,
                  volScalarField& weightField,
                  label myParticleType) const {}

  // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

  void resetVectorAverage(volVectorField& prev,
                          volVectorField& next,
                          bool single = false) const {
    if(!single) { prev == next; }
    next == dimensionedVector("zero", next.dimensions(), vector::zero);
  }

  void resetWeightFields() const {
    UsWeightField_ == dimensionedScalar("zero", UsWeightField_.dimensions(), 0.0);
  }

  void undoWeightFields(double** const& mask) const {
    for(int index = 0; index < particleCloud_.numberOfParticles(); index++) {
      // undo voidfraction cause by particle
      label cellI = particleCloud_.cellIDs()[index][0];
      UsWeightField_[cellI] -= particleCloud_.particleWeights()[index][0];
    }
  }

  tmp<volVectorField> UsInterp() const;

  virtual bool checkParticleType(label) const { return true; }

  // Access Functions
  volVectorField& UsPrev() const { return UsPrev_; }

  volVectorField& UsNext() const { return UsNext_; }

  volScalarField& UsWeightField() const { return UsWeightField_; }

};  // End of class averagingModel

}  // End of namespace Foam

#endif  // averagingModel_H
