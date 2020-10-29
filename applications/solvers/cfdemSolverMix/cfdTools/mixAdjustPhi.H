/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

InNamespace
    Foam

SourceFiles
    mixAdjustPhi.C

\*---------------------------------------------------------------------------*/

#ifndef mixAdjustPhi_H
#define mixAdjustPhi_H

#include "volFieldsFwd.H"
#include "surfaceFieldsFwd.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "inletOutletFvPatchFields.H"

namespace Foam {

//! @brief 在求解泊松方程的时候, 如果全部给定 Neumann 边界条件(即第二类边界条件), 还需要给定一个相容性条件, 这个附加条件就是要保证速度边界条件是守恒的
//        adjustPhi 和 pRef 针对的是压力方程为全Neuman边界的情况, 这会有两个问题:
//        (1) 方程系数矩阵 M 奇异, 解有无穷多, 相差一个常数, 所以需要 pRef
//        (2) 方程相容性问题, 需要 adjust Phi
//        比如方程为:
//          lapalacian(P) = f (in domain)
//          grad(P) & n = g (on boundary of domain, n 为边界单位外法向)
//        - pRef:
//          如果 p1 是上述方程的一个解，那么 p2 = p1 + const 也是方程的解，方程解不唯一，所以离散方程会有奇异性，无法求解
//        解决方法是固定一个点的压力为给定值，就是设置 pRef, 这个体现在主函数中的 pEqn.setReference(pRefCell, pRefValue)
//        - adjust Phi:
//          对于全Neumann条件的情形还需要满足相容性条件
bool mixAdjustPhi(surfaceScalarField& phi,
                  const volVectorField& U,
                  volScalarField& p) {

  if (!p.needReference()) {
    // 如果存在 Dirichlet 边界条件(第一类边界条件), 则 p.needReference() 为 false, 因此不需要满足附加条件
    return false;
  } else {

    scalar massIn = 0.0;
    scalar fixedMassOut = 0.0;
    scalar adjustableMassOut = 0.0;

    surfaceScalarField::Boundary& bphi = phi.boundaryFieldRef();

    forAll (bphi, patchi) {
      // 边界速度与通量
      const fvPatchVectorField& Up = U.boundaryField()[patchi];
      const fvsPatchScalarField& phip = bphi[patchi];
      
      if (!phip.coupled()) {
        if (Up.fixesValue() && !isA<inletOutletFvPatchVectorField>(Up)) {
          // 对于给定固定值且不是 inletOutletFvPatchVectorField 类型的速度边界, 计算总的入流量 massIn 以及总的固定出流量fixedMassOut
          forAll (phip, i) {
            if (phip[i] < 0.0) {
              massIn -= phip[i];
            }  else {
              fixedMassOut += phip[i];
            }
          }
        } else {
          // 对于其他速度边界, 计算总的入流量 massIn 和总的可调节出流量 adjustableMassOut
          forAll (phip, i) {
            if (phip[i] < 0.0) {
              massIn -= phip[i];
            } else {
              adjustableMassOut += phip[i];
            }
          }
        }
      }  // End of phip.coupled()
    }  // End of bphi loop

    // 计算当前计算域中总的通量(used for normalisation)
    scalar totalFlux = VSMALL + sum(mag(phi)).value();

    reduce(massIn, sumOp<scalar>());
    reduce(fixedMassOut, sumOp<scalar>());
    reduce(adjustableMassOut, sumOp<scalar>());

    scalar massCorr = 1.0;
    scalar magAdjustableMassOut = mag(adjustableMassOut);

    if (magAdjustableMassOut > VSMALL && magAdjustableMassOut / totalFlux > SMALL) {
      massCorr = (massIn - fixedMassOut) / adjustableMassOut;
    } else if (mag(fixedMassOut - massIn) / totalFlux > 1e-8) {
      // 当可调节流量很小，而流量的不平衡较大时报错
      FatalErrorInFunction
          << "Continuity error cannot be removed by adjusting the outflow.\nPlease check the velocity boundary conditions"
              " and/or run potentialFoam to initialise the outflow." << nl
          << "Total flux              : " << totalFlux << nl
          << "Specified mass inflow   : " << massIn << nl
          << "Specified mass outflow  : " << fixedMassOut << nl
          << "Adjustable mass outflow : " << adjustableMassOut << nl
          << exit(FatalError);
    }

    forAll (bphi, patchi) {
      const fvPatchVectorField& Up = U.boundaryField()[patchi];
      fvsPatchScalarField& phip = bphi[patchi];

      if (!phip.coupled()) {
        if (!Up.fixesValue() || isA<inletOutletFvPatchVectorField>(Up)) {
          forAll (phip, i) {
            if (phip[i] > 0.0) { phip[i] *= massCorr; }
          }
        }
      }  // End of phip.coupled()
    }  // End of bphi loop

    return mag(massIn) / totalFlux < SMALL
      && mag(fixedMassOut) / totalFlux < SMALL
      && mag(adjustableMassOut) / totalFlux < SMALL;
  }
}

}  // End of namespace Foam

#endif  // mixAdjustPhi_H
