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
\*---------------------------------------------------------------------------*/

#include "error.H"
#include "mixDense.H"
#include "addToRunTimeSelectionTable.H"
#include "voidFractionModel.H"

namespace Foam {

defineTypeNameAndDebug(mixDense, 0);

addToRunTimeSelectionTable(averagingModel, mixDense, dictionary);

// @brief Constructors
mixDense::mixDense(const dictionary& dict, cfdemCloud& sm): averagingModel(dict,sm) {}

// @brief Destructor
mixDense::~mixDense() {}

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
void mixDense::setScalarAverage(volScalarField& field,
                                double**& value,
                                double**& weight,
                                volScalarField& weightField,
                                double** const& mask,
                                double** const& weight2,
                                bool weightWithWeight2) const {
  FatalError << "mixDense::setScalarAverage() not implemented" << abort(FatalError);
}

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
void mixDense::setVectorAverage(volVectorField& field,
                                double**& value,
                                double**& weight,
                                volScalarField& weightField,
                                double** const& mask,
                                double** const& weight2,
                                bool weightWithWeight2) const {
  label cellI;
  vector valueVec;
  scalar weightP;

  for (int index = 0; index < particleCloud_.numberOfParticles(); index++) {

    //skip this particle if not correct type
    if (!checkParticleType(index)) { continue; }

    for (int subCell = 0; subCell < particleCloud_.cellsPerParticle()[index][0]; subCell++) {

      cellI = particleCloud_.cellIDs()[index][subCell];
      if (cellI >= 0) {
        for (int i = 0; i < 3; ++i) { valueVec[i] = value[index][i]; }
  
        if (weightWithWeight2) {
          weightP = weight[index][subCell] * weight2[index][subCell];
        } else {
          weightP = weight[index][subCell];
        }

        if(weightField[cellI] == 0) {  // first entry in this cell
          field[cellI] = valueVec;
          weightField[cellI] = weightP;
        } else {  // not first entry in this cell
          field[cellI] = (field[cellI] * weightField[cellI] + valueVec * weightP) / (weightField[cellI] + weightP);
          weightField[cellI] += weightP;
        }
      }  // End of cellI >= 0
    }  // End of subCell
  }  // End of index
  // correct cell values to patches
  field.correctBoundaryConditions();
}

/*!
 * \brief 设置局部平均矢量场
 * \param value              <[in] 用于局部平均化的颗粒(lagrange)变量
 * \param weight             <[in] 用于局部平均化的权重系数(lagrange)
 * \param dimensionRatios    <[in] 颗粒网格尺寸比
 * \param weightField        <[in, out] 权重系数平均化场
 * \param field              <[in, out] 需要被局部平均化场
 * \param mask
 */
void mixDense::setMixVectorAverage(double **& value,
                                   double **& weight,
                                   std::vector<double>& dimensionRatios,
                                   volVectorField& field,
                                   volScalarField& weightField,
                                   double **const& mask) const {
  label cellI;
  vector valueVec;
  scalar weightP;
  for (int index = 0; index < particleCloud_.numberOfParticles(); index++) {

    //skip this particle if not correct type
    if (!checkParticleType(index)) { continue; }

    if (!particleCloud_.useDynamicRefineMesh() && particleCloud_.checkFAndMParticle(dimensionRatios[index])) {

      for (int subcell = 0; subcell < particleCloud_.cellsPerParticle()[index][0]; ++subcell) {
        // 获取颗粒覆盖的第 subcell 个网格的 cellI
        cellI = particleCloud_.cellIDs()[index][subcell];

        if (cellI >= 0) {
          for (int i = 0; i < 3; ++i) { valueVec[i] = value[index][i]; }

          // 获取权重系数
          weightP = weight[index][subcell];

          if (weightField[cellI] == 0) {  // first entry in this cell
            field[cellI] = valueVec;
            weightField[cellI] = weightP;
          } else {  // not first entry in this cell
            field[cellI] = (field[cellI] * weightField[cellI] + valueVec * weightP) / (weightField[cellI] + weightP);
            weightField[cellI] += weightP;
          }
        }  // Cell found in domain
      }  // End of subcell
    }  // End of fine and middle particles
  }  // End of index

  // correct cell values to patches
  field.correctBoundaryConditions();
}

}  // End of namespace Foam

