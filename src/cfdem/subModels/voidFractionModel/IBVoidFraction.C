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
\*---------------------------------------------------------------------------*/

#include "subModels/voidFractionModel/IBVoidFraction.H"

namespace Foam {

cfdemDefineTypeName(IBVoidFraction)

cfdemCreateNewFunctionAdder(voidFractionModel, IBVoidFraction)

//! \brief Constructor
IBVoidFraction::IBVoidFraction(cfdemCloud& cloud)
    : voidFractionModel(cloud) {
}

//! \brief Destructor
IBVoidFraction::~IBVoidFraction() {}

/*!
 * \brief 设置单个颗粒的体积分数场
 * \param index 颗粒索引
 */
void IBVoidFraction::setVolumeFractionForSingleParticle(const int index) const {
  // 获取颗粒半径，中心坐标，颗粒中心所在网格编号
  double radius = cloud_.getRadius(index);
  Foam::vector centerPos = cloud_.getPosition(index);
}

/*!
 * \brief 计算距离系数，对任意一个网格, 如果网格中心 c 在颗粒内部, 但是它的某个角点 p
 *   不在颗粒内部, 则计算 c 与 p 的连线与颗粒表面的交点 i 到网格中心 c 的距离, 即
 *   求解 x 的二元一次方程
 *   (x * (vector_p - vector_c) - vector_particle) & 
 *   (x * (vector_p - vector_c) - vector_particle) == radius * radius
 *   等价于函数体中定义的: a*(x^2) - b*x + c = 0
 * \param index           <[in] 颗粒索引
 * \param pointInside     <[in] 网格中心
 * \param pointOutside    <[in] 网格角点
 */
double IBVoidFraction::segmentParticleIntersection(int index,
                                                   const Foam::vector& pointInside,
                                                   const Foam::vector& pointOutside) const {
  // 获取颗粒半径
  double radius = cloud_.getRadius(index);
  // 获取颗粒中心位置
  const Foam::vector centerPos = cloud_.getPosition(index);
  // 计算方程系数 a*(x^2) - b*x + c = 0
  double a = (pointOutside - pointInside) & (pointOutside - pointInside);
  double b = 2.0 * (pointOutside - pointInside) & (pointInside - centerPos);
  double c = ((pointInside - centerPos) & (pointInside - centerPos)) - radius * radius;
  double D = b * b - 4.0 * a * c;
  double lambda_1 = 0.0;
  double lambda_2 = 0.0;
  double lambda = 0.0;
  double eps = 1.0e-15;
  if (D >= 0.0) {
    // 方程有实数解, 则说明连线与球面一定有两个交点, 分别用 lambda_1 和 lambda_2 表示方程的两个实数解
    lambda_1 = (-1.0 * b + sqrt(D)) / (2.0 * a);
    lambda_2 = (-1.0 * b - sqrt(D)) / (2.0 * a);
    if (lambda_1 >= -eps && lambda_1 <= 1.0 + eps) {
      // 如果 lambda_1 在 [0, 1] 之间, 则说明这个交点在网格中心与网格角点的连线内, 则返回这个解
      lambda = lambda_1;
    } else if (lambda_2 >= -eps && lambda_2 <= 1.0 + eps) {
      // 如果 lambda_2 在 [0, 1] 之间, 则说明这个交点在网格中心与网格角点的连线内, 则返回这个解
      lambda = lambda_2;
    }
  }
  // 确保 lambda 在 [0, 1] 区间内, 而不是在 [-eps, 1 + eps] 区间内
  if (lambda > 1.0) {
    return 1;
  }
  if (lambda < 0.0) {
    return 0.0;
  }
  return lambda;
}

} // namespace Foam
