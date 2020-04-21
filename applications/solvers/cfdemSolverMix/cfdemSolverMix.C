/*---------------------------------------------------------------------------*\
    CFDEMcoupling - Open Source CFD-DEM coupling

    CFDEMcoupling is part of the CFDEMproject
    www.cfdem.com
                                Christoph Goniva, christoph.goniva@cfdem.com
                                Copyright (C) 1991-2009 OpenCFD Ltd.
                                Copyright (C) 2009-2012 JKU, Linz
                                Copyright (C) 2012-     DCS Computing GmbH,Linz
-------------------------------------------------------------------------------
License
    This file is part of CFDEMcoupling.

    CFDEMcoupling is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    CFDEMcoupling is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with CFDEMcoupling.  If not, see <http://www.gnu.org/licenses/>.

Application
    cfdemSolverMix

Description
    cfdemSolverMix
\*---------------------------------------------------------------------------*/
#ifndef __NONETYPE__
#define __NONETYPE__ 1
#endif

#include "fvCFD.H"
#include "singlePhaseTransportModel.H"

#include "OFversion.H"
#if defined(version30)
  #include "turbulentTransportModel.H"
  #include "pisoControl.H"
#else
  #include "turbulenceModel.H"
#endif

#if defined(versionv1606plus) || defined(version40)
  #include "fvOptions.H"
#else
  #include "fvIOoptionList.H"
#endif

#include "fixedFluxPressureFvPatchScalarField.H"
#include "../cfdemSolverMix/cfdemCloudMix.H"

#if defined(anisotropicRotation)
  #include "cfdemCloudRotation.H"
#endif

#if defined(superquadrics_flag)
  #include "cfdemCloudRotationSuperquadric.H"
#endif

#include "implicitCouple.H"
#include "clockModel.H"
#include "smoothingModel.H"
#include "forceModel.H"

#include "./cfdTools/mixAdjustPhi.H"
#include "dynamicFvMesh.H"

int main(int argc, char *argv[]) {
  #include "setRootCase.H"
  #include "createTime.H"

  // #include "createMesh.H"
  #include "createDynamicFvMesh.H"

#if defined(version30)
  pisoControl piso(mesh);
  #include "createTimeControls.H"
#endif

  // 创建所有场变量
  #include "createFields.H"

  #include "createFvOptions.H"
  #include "initContinuityErrs.H"

  #include "readGravitationalAcceleration.H"
  // 在 cfdemSolverMix 求解器中必须使用 impCoupleModel
  // 在 cfdemSolverMix 求解器中必须使用 mix force model
  #include "./mixCheckGlobal.H"

  // create cfdemCloud
#if defined(anisotropicRotation)
  cfdemCloudRotation particleCloud(mesh);
#elif defined(superquadrics_flag)
  cfdemCloudRotationSuperquadric particleCloud(mesh);
#else
  cfdemCloudMix particleCloud(mesh);
#endif

  #include "./mixCheckModelType.H"
  word modelType = particleCloud.modelType();

  // * * * * * * * * * * * * * * * * * Starting time loop * * * * * * * * * * * * * * * * * //
  Info << "\nStarting time loop\n" << endl;

  while(runTime.loop()) {
  
    Info << "Time = " << runTime.timeName() << nl << endl;

    // 设置动态加密网格
    interFace = mag(mesh.lookupObject<volScalarField>("volumefractionNext"));
    particleCloud.setMeshHasUpdatedFlag(mesh.update());

#if defined(version30)
    #include "readTimeControls.H"
    #include "CourantNo.H"
    #include "setDeltaT.H"
#else
    #include "readPISOControls.H"
    #include "CourantNo.H"
#endif

    bool hasEvolved = particleCloud.evolve(voidfraction, volumefraction, Us, U, interFace);

    if (hasEvolved) {
      // 光滑隐式力场
      particleCloud.smoothingM().smoothenAbsolutField(particleCloud.forceM(0).impParticleForces());
    }

    // 获取单位体积动量交换场
    // 在 mixCheckImCoupleM.H 中限制了 momCoupleModel 只能为隐式交换, 所以在 getProperty 函数中只能是: "implicitCouple_index" or "mixImplicitCouple_index"
    Ksl = particleCloud.momCoupleM(particleCloud.registryM().getProperty("mixImplicitCouple_index")).impMomSource();
    Ksl.correctBoundaryConditions();

    // #include "forceCheckIm.H"
    {
      Info << "cfdemSolverMix: forceCheckIm..." << endl;
      // 单位体积动量交换系数与局部平均颗粒速度都是当前时间步的值, 流体速度则是上个时间步的值, 定义单位体积隐式力 fImp = Ksl * (Us - U)
      volVectorField fImp(Ksl * (Us - U));
      // 计算隐式力的值: forAll(fImp, cellI) { fImp[cellI] *= particleCloud.mesh().V()[cellI]; }
      particleCloud.scaleWithVcell(fImp);
      dimensionedVector fImpTotal = gSum(fImp);
      Info << "TotalForceImp:  " << fImpTotal.value() << endl;
      Info << "Warning, these values are based on latest Ksl and Us but prev. iteration U!" << endl;
      Info << "cfdemSolverMix: forceCheckIm - done\n" << endl;
    }
    #include "solverDebugInfo.H"

    // 如果 modelType 为 none, 则说明当前计算忽略所有的 fine particles
    #include "./mixCheckModelNone.H"

#if __NONETYPE__
    surfaceScalarField voidfractionf = fvc::interpolate(voidfraction);
#else
    surfaceScalarField voidfractionf = fvc::interpolate(voidfraction);
    phi = voidfractionf * phiByVoidfraction;
#endif

    if (particleCloud.solveFlow()) {

#if __NONETYPE__
      fvVectorMatrix UEqn
      (
          fvm::ddt(U)
        + fvm::div(phi, U)
        + turbulence->divDevReff(U)
        ==
          fvOptions(U)
      );
#else
      // 动量预测
      // 注意: 在动量方程中, modelType 为 "B" or "Bfull" 的时候, 应力项中不需要乘以空隙率, 而当 modelType 为 "A" 时, 应力项中需要乘以空隙率
      fvVectorMatrix UEqn
      (
          fvm::ddt(voidfraction, U)
        - fvm::Sp(fvc::ddt(voidfraction), U)
        + fvm::div(phi, U)
        - fvm::Sp(fvc::div(phi), U)
        // + particleCloud.divVoidfractionTau(U, voidfraction)
        - fvm::laplacian(particleCloud.voidfractionNuEff(voidfraction), U)
        - fvc::div(particleCloud.voidfractionNuEff(voidfraction) * dev2(fvc::grad(U)().T()))
        // + turbulence->divDevReff(U)
        ==
        - fvm::Sp(Ksl / rho, U)
        + fvOptions(U)
      );
#endif  // __NONETYPE__

      UEqn.relax();
      fvOptions.constrain(UEqn);

#if defined(version30)
      if (piso.momentumPredictor())
#else
      if (momentumPredictor)
#endif
      {
        if (modelType == "B" || modelType == "Bfull") {
          // 在动量方程中, modelType 为 "B" or "Bfull" 的时候, 压力项中不需要乘以空隙率
          solve(UEqn == - fvc::grad(p) + Ksl / rho * Us);
        } else if (modelType == "A") {
          // 在动量方程中, modelType 为 "A" 时, 压力项中需要乘以空隙率
          solve(UEqn == - voidfraction * fvc::grad(p) + Ksl / rho * Us);
        } else {  // modelType == "none"
#if __NONETYPE__
          solve(UEqn == -fvc::grad(p));
#else
          FatalError << "cfdemSolverMix.C: " << __LINE__ << ": Not implement for modelType = None" << abort(FatalError);
#endif  // __NONETYPE__
        }
        fvOptions.correct(U);
      }  // End of momentumPredictor

      // PISO loop
#if defined(version30)
      while (piso.correct())
#else
      for (int corr = 0; corr < nCorr; corr++)
#endif
      {
        volScalarField rUA = 1.0 / UEqn.A();
        volVectorField HbyA = rUA * UEqn.H();

        if (modelType != "none") {  // "A" or "B" or  "Bfull"

          // - 定义 HbyA 的通量场(不包含 Ksl / rho * Us 的通量)
          surfaceScalarField phiHbyA(
            "phiHbyA",
#ifdef version23
            (fvc::interpolate(HbyA) & mesh.Sf()) + fvc::interpolate(rUA) * fvc::ddtCorr(HbyA, phiByVoidfraction)
#else
            (fvc::interpolate(HbyA) & mesh.Sf()) + fvc::ddtPhiCorr(rUA, HbyA, phiByVoidfraction);
#endif
          );

          // - 定义 Ksl / rho * Us 的通量
          surfaceScalarField phiUs(fvc::interpolate(Us) & mesh.Sf());

          // - 将Us通量与HbyA通量相加，得到压力修正方程中的总通量
          phiHbyA += fvc::interpolate(rUA) * fvc::interpolate(Ksl / rho) * phiUs;

          // - 修正通量 phiHbyA 以满足泊松方程相容性条件(当压力边界全部为 Neumann 条件)
          // mixAdjustPhi(phiHbyA, U, p); ???????????????

          // - 通过迭代求解压力泊松方程, 并对方程中的拉普拉斯项进行非正交修正
          // 注意: 由于压力泊松方程中存在压力的拉普拉斯项, 当网格非正交的时候会出现显式源项, 显式源项会在正交修正迭代后进行更新
          // 修正次数可由 controlDict 中的 nNonOrthogonalCorrectors 来指定
#if defined(version30)
          while (piso.correctNonOrthogonal())
#else 
          for (int nonOrth = 0; nonOrth <= nNonOrthCorr; ++nonOrth)
#endif
          {
            volScalarField rUAvoidfraction("(voidfraction / A(U))", rUA * voidfraction);
            if (modelType == "A") {
              rUAvoidfraction = volScalarField("voidfraction^2 / A(U)", rUA * voidfraction * voidfraction);
            }

            fvScalarMatrix pEqn
            (
                fvm::laplacian(rUAvoidfraction, p)
              ==
                fvc::div(voidfractionf * phiHbyA) + particleCloud.ddtVoidfraction()
            );
            pEqn.setReference(pRefCell, pRefValue);

#if defined(version30)

            pEqn.solve(mesh.solver(p.select(piso.finalInnerIter())));
            if (piso.finalNonOrthogonalIter()) {
              phiByVoidfraction = phiHbyA - pEqn.flux() / voidfractionf;
            }

#else

            if (corr == nCorr-1 && nonOrth == nNonOrthCorr) {
  #if defined(versionExt32)
              pEqn.solve(mesh.solutionDict().solver("pFinal"));
  #else
              pEqn.solve(mesh.solver("pFinal"));
  #endif
            } else {
              pEqn.solve();
            }
            if(nonOrth == nNonOrthCorr) {
              phiByVoidfraction = phiHbyA - pEqn.flux() / voidfractionf;
            }

#endif
          }  // End of non-orthogonal corrector loop

          phi = voidfractionf * phiByVoidfraction;
          #include "continuityErrorPhiPU.H"

          // 更新速度场
          if (modelType == "A") {
            U = HbyA - rUA * (voidfraction * fvc::grad(p) - Ksl / rho * Us);
          } else if (modelType == "B" || modelType == "Bfull") {
            U = HbyA - rUA * (fvc::grad(p) - Ksl / rho * Us);
          }
          U.correctBoundaryConditions();
          fvOptions.correct(U);

        } else {  // "none"
#if __NONETYPE__
          surfaceScalarField phiHbyA(
            "phiHbyA",
#ifdef version23
            (fvc::interpolate(HbyA) & mesh.Sf()) + fvc::interpolate(rUA) * fvc::ddtCorr(HbyA, phi)
#else
            (fvc::interpolate(HbyA) & mesh.Sf()) + fvc::ddtPhiCorr(rUA, HbyA, phi);
#endif
          );
          adjustPhi(phiHbyA, U, p);

#if defined(version30)
          while (piso.correctNonOrthogonal())
#else 
          for (int nonOrth = 0; nonOrth <= nNonOrthCorr; ++nonOrth)
#endif
          {
            fvScalarMatrix pEqn
            (
              fvm::laplacian(rUA, p) == fvc::div(phiHbyA) + particleCloud.ddtVoidfraction()
            );
            pEqn.setReference(pRefCell, pRefValue);

#if defined(version30)
            pEqn.solve(mesh.solver(p.select(piso.finalInnerIter())));
            if (piso.finalNonOrthogonalIter()) {
              phi = phiHbyA - pEqn.flux();
            }
#else
            if (corr == nCorr-1 && nonOrth == nNonOrthCorr) {
  #if defined(versionExt32)
              pEqn.solve(mesh.solutionDict().solver("pFinal"));
  #else
              pEqn.solve(mesh.solver("pFinal"));
  #endif
            } else {
              pEqn.solve();
            }
            if (nonOrth == nNonOrthCorr) {
              phi = phiHbyA - pEqn.flux();
            }
#endif
          }
          #include "continuityErrs.H"
          U = HbyA - rUA * fvc::grad(p);
          U.correctBoundaryConditions();
#else
          FatalError << "cfdemSolverMix.C: " << __LINE__ << ": Not implement for modelType = None" << abort(FatalError);
#endif  // __NONETYPE__
        }

      }  // End of PISO loop

      laminarTransport.correct();
      turbulence->correct();

#if __NONETYPE__
      volScalarField volumefractionNext = mesh.lookupObject<volScalarField>("volumefractionNext");
      Info << "cfdemSolverMix: calcVelocityCorrection..." << endl;
      particleCloud.calcVelocityCorrection(p, U, phiIB, volumefractionNext);
      Info << "cfdemSolverMix: calcVelocityCorrection - done" << endl;
  #if defined(version22)
      fvOptions.correct(U);
  #endif
#endif  // __NONETYPE__

    // End of solveFlow()
    } else {
      Info << "cfdemSolverMix.C: " << __LINE__ << ": Skipping flow solution." << endl;
    }

    runTime.write();
    Info << "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
      << "  ClockTime = " << runTime.elapsedClockTime() << " s" << nl << endl;

  }  // End of while(runTime.loop())

  Info << "cfdemSolverMix.C: End\n" << endl;
  return 0;
}
