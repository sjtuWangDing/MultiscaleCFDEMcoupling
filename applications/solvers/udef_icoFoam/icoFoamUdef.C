#include "fvCFD.H"
#include "pisoControl.H"

int main(int argc, char *argv[]) {
  #include "setRootCase.H"
  #include "createTime.H"
  #include "createMesh.H"

  pisoControl piso(mesh);

  #include "createFields.H"
  #include "initContinuityErrs.H"

  Info << "\nStarting time loop\n" << endl;
  while (runTime.loop()) {
    Info << "Time = " << runTime.timeName() << nl << endl;
    #include "CourantNo.H"
    // Momentum predictor
    fvVectorMatrix UEqn(
        fvm::ddt(U)
      + fvm::div(phi, U)
      - fvm::laplacian(nu, U)
    );

    if (piso.momentumPredictor()) {
      solve(UEqn == -fvc::grad(p));
    }

    // PISO loop
    while (piso.correct()) {

      // UEqn.A() 为方程(13)中 Ap
      volScalarField rAU(1.0 / UEqn.A());
      volVectorField HbyA(constrainHbyA(rAU * UEqn.H(), U, p));  // 也可以写为 HbyA = rAU * UEqn.H()

      // phiHbyA 为方程(25)中HybA的通量, 等于主通量 + 修正通量
      surfaceScalarField phiHbyA
      (
          "phiHbyA",
          fvc::flux(HbyA) + fvc::interpolate(rAU) * fvc::ddtCorr(U, phi)
          // (fvc::interpolate(HbyA) & mesh.Sf()) + fvc::ddtPhiCorr(rAU, U, phi)
      );
      adjustPhi(phiHbyA, U, p);

      // Update the pressure BCs to ensure flux consistency
      constrainPressure(p, U, phiHbyA, rAU);

      // Non-orthogonal pressure corrector loop
      while (piso.correctNonOrthogonal()) {
        // Pressure corrector
        fvScalarMatrix pEqn
        (
          fvm::laplacian(rAU, p) == fvc::div(phiHbyA)
        );
        pEqn.setReference(pRefCell, pRefValue);

        pEqn.solve(mesh.solver(p.select(piso.finalInnerIter())));

        if (piso.finalNonOrthogonalIter()) {
          // finalNonOrthogonalIter() return corrNonOrtho_ == nNonOrthCorr_ + 1;
          phi = phiHbyA - pEqn.flux();
        }
      }  // End of pressure corrector loop

      // #include "continuityErrs.H"
      {
        volScalarField contErr(fvc::div(phi));

        scalar sumLocalContErr = runTime.deltaTValue() * mag(contErr)().weightedAverage(mesh.V()).value();

        scalar globalContErr = runTime.deltaTValue() * contErr.weightedAverage(mesh.V()).value();
        cumulativeContErr += globalContErr;

        Info << "time step continuity errors : sum local = " << sumLocalContErr
             << ", global = " << globalContErr
             << ", cumulative = " << cumulativeContErr << endl;
      }
      U = HbyA - rAU * fvc::grad(p);
      U.correctBoundaryConditions();
    }  // End of piso loop

    runTime.write();
    Info << "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
         << "  ClockTime = " << runTime.elapsedClockTime() << " s"
         << nl << endl;
  }  // End of runTime.loop()

  Info<< "End\n" << endl;
  return 0;
}