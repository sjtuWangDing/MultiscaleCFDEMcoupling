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
#include "cfdemCloud.H"

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

int main(int argc, char *argv[]) {
  #include "setRootCase.H"
  #include "createTime.H"
  #include "createMesh.H"

#if defined(version30)
  pisoControl piso(mesh);
  #include "createTimeControls.H"
#endif

  #include "createFields.H"
  #include "createFvOptions.H"
  #include "initContinuityErrs.H"

  #include "readGravitationalAcceleration.H"
  #include "checkImCoupleM.H"
  // Info << "\nReading g" << endl;
  // uniformDimensionedVectorField g
  // (
  //   IOobject
  //   (
  //     "g",
  //     runTime.constant(),
  //     mesh,
  //     IOobject::MUST_READ,
  //     IOobject::NO_WRITE
  //   )
  // );

  // create cfdemCloud
#if defined(anisotropicRotation)
  cfdemCloudRotation particleCloud(mesh);
#elif defined(superquadrics_flag)
  cfdemCloudRotationSuperquadric particleCloud(mesh);
#else
  cfdemCloud particleCloud(mesh);
#endif
  // check model "A" "B" "Bfull"
  #include "checkModelType.H"

  Info << "\nStarting time loop\n" << endl;
  while (runTime.loop()) {
    Info << "Time = " << runTime.timeName() << nl << endl;

#if defined(version30)
    #include "readTimeControls.H"
    #include "CourantNo.H"
    #include "setDeltaT.H"
#else
    #include "readPISOControls.H"
    #include "CourantNo.H"
#endif

  bool hasEvolved = particleCloud.evolve(voidfraction, Us, U);

    if (hasEvolved) {
      particleCloud.smoothingM().smoothenAbsolutField(particleCloud.forceM(0).impParticleForces());
    }
    Ksl = particleCloud.momCoupleM(particleCloud.registryM().getProperty("implicitCouple_index")).impMomSource();
    Ksl.correctBoundaryConditions();
    surfaceScalarField voidfractionf = fvc::interpolate(voidfraction);
    phi = voidfractionf * phiByVoidfraction;

    // #include "forceCheckIm.H"
    {
      Info << "\nSolver level total Eulerian momentum exchange:"<< endl;
      // calculate total implicit force
      volVectorField fImp(Ksl * (Us - U));
      particleCloud.scaleWithVcell(fImp);
      dimensionedVector fImpTotal = gSum(fImp);
      Info << "  TotalForceImp:  " << fImpTotal.value() << endl;
      Info << "  Warning, these values are based on latest Ksl and Us but prev. iteration U!\n" << endl;
    }
    #include "solverDebugInfo.H"

    if(particleCloud.solveFlow()) { // Pressure-velocity PISO corrector
      // Momentum predictor
      fvVectorMatrix UEqn
      (
          fvm::ddt(voidfraction, U) - fvm::Sp(fvc::ddt(voidfraction), U)
        + fvm::div(phi, U) - fvm::Sp(fvc::div(phi), U)
        // + turbulence->divDevReff(U)
        + particleCloud.divVoidfractionTau(U, voidfraction)
        ==
        - fvm::Sp(Ksl/rho,U)
        + fvOptions(U)
      );
    }
  }  // End of runTime.loop().
}