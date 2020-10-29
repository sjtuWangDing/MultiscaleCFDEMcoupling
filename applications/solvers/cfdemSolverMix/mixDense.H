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

    averaging model for dense regime

Class
    mixDense

SourceFiles
    mixDense.C
\*---------------------------------------------------------------------------*/
#ifndef mixDense_H
#define mixDense_H

#include "averagingModel.H"

namespace Foam {

class mixDense: public averagingModel {

public:

  //! @brief Runtime type information
  TypeName("mixDense");

  //! @brief Constructors
  mixDense(const dictionary& dict, cfdemCloud& sm);

  //! @brief Destructor
  ~mixDense();

  /*!
   * \brief 设置局部平均标量场
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
                                bool weightWithWeight2 = false) const;

  /*!
   * \brief 设置局部平均矢量场
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
                                bool weightWithWeight2 = false) const;

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
                                   double **const& mask) const;

};  // class mixDense

}  // namespace Foam

#endif  // mixDense_H
