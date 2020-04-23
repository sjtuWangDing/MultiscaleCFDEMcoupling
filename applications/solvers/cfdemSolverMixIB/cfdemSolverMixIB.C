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

// #include "./cfdTools/mixAdjustPhi.H"
#include "dynamicFvMesh.H"

int main(int argc, char *argv[]) {

  #include "setRootCase.H"
  #include "createTime.H"
  #include "createDynamicFvMesh.H"

#if defined(version30)
  pisoControl piso(mesh);
  #include "createTimeControls.H"
#endif

  #include "./createFields.H"
  #include "initContinuityErrs.H"
#if defined(version22)
  #include "createFvOptions.H"
#endif

  #include "readGravitationalAcceleration.H"
  // create cfdemCloud
#if defined(anisotropicRotation)
  cfdemCloudRotation particleCloud(mesh);
#elif defined(superquadrics_flag)
  cfdemCloudRotationSuperquadric particleCloud(mesh);
#else
  cfdemCloudMix particleCloud(mesh);
#endif

  word modelType = particleCloud.modelType();
  if (modelType != "none") {
    FatalError << "cfdemSolverMixIB.C: modelType is " << modelType << " ,not none" << abort(FatalError);
  }

  // * * * * * * * * * * * * * * * * * Starting time loop * * * * * * * * * * * * * * * * * //
  Info << "\nStarting time loop\n" << endl;
  while (runTime.loop()) {
    Info << "Time = " << runTime.timeName() << nl << endl;

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
    if (particleCloud.solveFlow()) {

      fvVectorMatrix UEqn
      (
          fvm::ddt(U)
        + fvm::div(phi, U)
        + turbulence->divDevReff(U)
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

  // --- PISO loop
#if defined(version30)
      while (piso.correct())
#else
      for (int corr = 0; corr < nCorr; corr++)
#endif
      {
        volScalarField rUA = 1.0 / UEqn.A();
        surfaceScalarField rUAf(fvc::interpolate(rUA));
        U = rUA * UEqn.H();
#ifdef version23
        phi = (fvc::interpolate(U) & mesh.Sf()) + rUAf*fvc::ddtCorr(U, phi);
#else
        phi = (fvc::interpolate(U) & mesh.Sf()) + fvc::ddtPhiCorr(rUA, U, phi);
#endif
        adjustPhi(phi, U, p);

#if defined(version22)
        fvOptions.relativeFlux(phi);
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
            fvm::laplacian(rUA, p) == fvc::div(phi) + particleCloud.ddtVoidfraction()
          );
          pEqn.setReference(pRefCell, pRefValue);

#if defined(version30)
          pEqn.solve(mesh.solver(p.select(piso.finalInnerIter())));
          if (piso.finalNonOrthogonalIter()) {
            phi -= pEqn.flux();
          }
#else
          if (corr == nCorr-1 && nonOrth == nNonOrthCorr ) {
  #if defined(versionExt32)
            pEqn.solve(mesh.solutionDict().solver("pFinal"));
  #else
            pEqn.solve(mesh.solver("pFinal"));
  #endif
          } else {
            pEqn.solve();
          }

          if (nonOrth == nNonOrthCorr) {
            phi -= pEqn.flux();
          }
#endif
        }
        #include "continuityErrs.H"
        U -= rUA * fvc::grad(p);
        U.correctBoundaryConditions();
      }  // End of piso loop
    }  // End of particleCloud.solveFlow

    laminarTransport.correct();
    turbulence->correct();

    volScalarField volumefractionNext = mesh.lookupObject<volScalarField>("volumefractionNext");
    Info << "particleCloud.calcVelocityCorrection()..." << endl;
    particleCloud.calcVelocityCorrection(p, U, phiIB, volumefractionNext);
    Info << "particleCloud.calcVelocityCorrection() - done" << endl;

#if defined(version22)
    fvOptions.correct(U);
#endif

    runTime.write();

    Info << "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
      << "  ClockTime = " << runTime.elapsedClockTime() << " s" << nl << endl;
  }
  Info << "End\n" << endl;
  return 0;
}
