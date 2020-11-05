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
  #include "createTimeControls.H"
  #include "readGravitationalAcceleration.H"
  pisoControl piso(mesh);

  #include "./createFields.H"
  #include "initContinuityErrs.H"
#if defined(version22)
  #include "createFvOptions.H"
#endif

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

    #include "readTimeControls.H"
    #include "CourantNo.H"
    #include "setDeltaT.H"

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

      if (piso.momentumPredictor()) {
        solve(UEqn == -fvc::grad(p));
      }

      // PISO loop
      while (piso.correct()) {
        // 定义 HbyA，计算并修正 HbyA 的通量 phiHbyA
        volScalarField rAU(1.0 / UEqn.A());
        volVectorField HbyA(constrainHbyA(rAU * UEqn.H(), U, p));
        surfaceScalarField phiHbyA(
          "phiHbyA",
          fvc::flux(HbyA) + fvc::interpolate(rAU) * fvc::ddtCorr(U, phi)
        );

        // 在 adjustPhi 函数中，第二个参数必须使用 U，而不能使用 HbyA，因为在 adjustPhi 函数中，需要通过 U 获取 boundaryField
        adjustPhi(phiHbyA, U, p);

        // Update the pressure BCs to ensure flux consistency
        constrainPressure(p, U, phiHbyA, rAU);

        // Non-orthogonal pressure corrector loop
        while (piso.correctNonOrthogonal()) {
          // Pressure corrector
          fvScalarMatrix pEqn (
            fvm::laplacian(rAU, p) == fvc::div(phiHbyA) + particleCloud.ddtVoidfraction()
          );
          pEqn.setReference(pRefCell, pRefValue);
          pEqn.solve(mesh.solver(p.select(piso.finalInnerIter())));

          if (piso.finalNonOrthogonalIter()) {
            phi = phiHbyA - pEqn.flux();
          }
        }
        #include "continuityErrs.H"
        U = HbyA - rAU * fvc::grad(p);
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
