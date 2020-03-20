#include "fvCFD.H"
#include "singlePhaseTransportModel.H"
#include "OFversion.H"

#if defined(version30)
  #include "turbulentTransportModel.H"
  #include "pisoControl.H"
#else
  #include "turbulenceModel.H"
#endif

#include "./udef_cfdemCloudIB.H"

#if defined(superquadrics_flag)
  #include "cfdemCloudIBSuperquadric.H"
#endif

#include "implicitCouple.H"
#include "averagingModel.H"
#include "voidFractionModel.H"
#include "dynamicFvMesh.H"
#include "cellSet.H"

#if defined(version22)
  #include "meshToMeshNew.H"
  #include "fvIOoptionList.H"
#endif

int main(int argc, char *argv[]) {
  #include "setRootCase.H"          // 设置算例的根目录
  #include "createTime.H"           //创建时间对象
  #include "createDynamicFvMesh.H"  //创建网格对象

  #if defined(version30)
    pisoControl piso(mesh);
    #include "createTimeControls.H"
  #endif

  #include "createFields.H"
  #include "initContinuityErrs.H"

  #if defined(version22)
    #include "createFvOptions.H"
  #endif

  // create cfdemCloud
  #include "readGravitationalAcceleration.H"
  #if defined(superquadrics_flag)
    cfdemCloudIBSuperquadric particleCloud(mesh);
  #else
    udef_cfdemCloudIB particleCloud(mesh);
  #endif

  Info << "\nStarting time loop\n" << endl;

  while (runTime.loop()) {
    Info << "Time = " << runTime.timeName() << nl << endl;
    // dynamic grid
    interFace = mag(mesh.lookupObject<volScalarField>("voidfractionNext"));
    particleCloud.setMeshHasUpdatedFlag(mesh.update());

    #if defined(version30)
      #include "readTimeControls.H"
      #include "CourantNo.H"
      #include "setDeltaT.H"
    #else
      #include "readPISOControls.H"
      #include "CourantNo.H"
    #endif

    Info << "- evolve()" << endl;
    particleCloud.evolve(voidfraction, interFace);

    // Pressure-velocity PISO corrector
    if (particleCloud.solveFlow()) {
      // Momentum predictor
      fvVectorMatrix UEqn
      (
          fvm::ddt(U)
        + fvm::div(phi, U)
        // + turbulence->divDevReff(U)
        #if defined(version22)
        ==
          fvOptions(U)
        #endif
      );
      UEqn.relax();
#if defined(version22)
      fvOptions.constrain(UEqn);
#endif

#if defined(version30)
      if (piso.momentumPredictor())
#else
      if (momentumPredictor)
#endif
      {
        solve(UEqn == -fvc::grad(p));
      }

      // - piso loop
#if defined(version30)
      while (piso.correct())
#else
      for (int corr = 0; corr < nCorr; corr++)
#endif
      {
        volScalarField rAU = 1.0 / UEqn.A();
        volVectorField HbyA = rAU * UEqn.H();
        surfaceScalarField phiHbyA
        (
          "phiHbyA",  // ??????????????
#ifdef version23
          (fvc::interpolate(HbyA) & mesh.Sf()) + fvc::interpolate(rAU) * fvc::ddtCorr(U, phi)
#else
          (fvc::interpolate(HbyA) & mesh.Sf()) + fvc::ddtPhiCorr(rAU, U, phi);
#endif
        );
        adjustPhi(phiHbyA, U, p);
        #if defined(version22)
        fvOptions.relativeFlux(phiHbyA);
        #endif

        // Non-orthogonal pressure corrector loop
#if defined(version30)
        while (piso.correctNonOrthogonal())
#else
        for (int nonOrth = 0; nonOrth <= nNonOrthCorr; nonOrth++)
#endif
        {
          // Pressure corrector
          fvScalarMatrix pEqn
          (
            fvm::laplacian(rAU, p) == fvc::div(phiHbyA)
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
        }  // End of non-orthogonal pressure corrector loop
        U = HbyA - rAU * fvc::grad(p);
        U.correctBoundaryConditions();
      }  // End of piso loop
    }  // End of solveFlow
    laminarTransport.correct();
    turbulence->correct();

    Info << "particleCloud.calcVelocityCorrection() " << endl;
    volScalarField voidfractionNext = mesh.lookupObject<volScalarField>("voidfractionNext");
    particleCloud.calcVelocityCorrection(p, U, phiIB, voidfractionNext);

#if defined(version22)
    fvOptions.correct(U);
#endif

    runTime.write();
    Info << "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
      << "  ClockTime = " << runTime.elapsedClockTime() << " s"
      << nl << endl;
  }  // End of runTime loop
  Info << "End\n" << endl;
  return 0;
}