Info << "Reading field p\n" << endl;
volScalarField p
(
  IOobject
  (
    "p",
    runTime.timeName(),
    mesh,
    IOobject::MUST_READ,
    IOobject::AUTO_WRITE
  ),
  mesh
);

Info << "Reading physical velocity field U" << endl;
Info << "Note: only if voidfraction at boundary is 1, U is superficial velocity!!!\n" << endl;
volVectorField U
(
  IOobject
  (
    "U",
    runTime.timeName(),
    mesh,
    IOobject::MUST_READ,
    IOobject::AUTO_WRITE
  ),
  mesh
);

Info << "Reading particle velocity field Us\n" << endl;
volVectorField Us
(
  IOobject
  (
    "Us",
    runTime.timeName(),
    mesh,
    IOobject::READ_IF_PRESENT,
    IOobject::AUTO_WRITE
  ),
  mesh,
  dimensionedVector("zero", dimensionSet(0, 1, -1, 0, 0), vector(0, 0, 0))
);

Info << "\nReading momentum exchange field Ksl\n" << endl;
volScalarField Ksl
(
  IOobject
  (
    "Ksl",
    runTime.timeName(),
    mesh,
    IOobject::READ_IF_PRESENT,
    IOobject::AUTO_WRITE
  ),
  mesh,
  dimensionedScalar("zero", dimensionSet(1, -3, -1, 0, 0), 0)
);

Info << "\nReading voidfraction field voidfraction = (Vgas/Vparticle)\n" << endl;
volScalarField voidfraction
(
  IOobject
  (
    "voidfraction",
    runTime.timeName(),
    mesh,
    IOobject::READ_IF_PRESENT,
    IOobject::AUTO_WRITE
  ),
  mesh,
  dimensionedScalar("zero", dimensionSet(0, 0, 0, 0, 0), 1)
);

Info << "\nReading volumefraction field volumefraction = (Vgas/Vparticle)\n" << endl;
volScalarField volumefraction
(
  IOobject
  (
    "volumefraction",
    runTime.timeName(),
    mesh,
    IOobject::READ_IF_PRESENT,
    IOobject::AUTO_WRITE
  ),
  mesh,
  dimensionedScalar("zero", dimensionSet(0, 0, 0, 0, 0), 1)
);

Info << "\nReading density field rho\n" << endl;
volScalarField rho
(
  IOobject
  (
    "rho",
    runTime.timeName(),
    mesh,
    IOobject::MUST_READ,
    IOobject::AUTO_WRITE
  ),
  mesh
);

#ifndef createPhi_H
#define createPhi_H
Info << "Reading/calculating face flux field phi\n" << endl;
surfaceScalarField phi
(
  IOobject
  (
    "phi",
    runTime.timeName(),
    mesh,
    IOobject::READ_IF_PRESENT,
    IOobject::AUTO_WRITE
  ),
  linearInterpolate(U * voidfraction) & mesh.Sf()
);
#endif

Info << "Generating interstitial face flux field phiByVoidfraction\n" << endl;
surfaceScalarField phiByVoidfraction
(
  IOobject
  (
    "phiByVoidfraction",
    runTime.timeName(),
    mesh,
    IOobject::NO_READ,
    IOobject::NO_WRITE
  ),
  linearInterpolate(U) & mesh.Sf()
);

label pRefCell = 0;
scalar pRefValue = 0.0;
setRefCell(p, mesh.solutionDict().subDict("PISO"), pRefCell, pRefValue);

singlePhaseTransportModel laminarTransport(U, phi);

autoPtr<incompressible::turbulenceModel> turbulence(
  incompressible::turbulenceModel::New(U, phi, laminarTransport)
);

Info << "Reading field phiIB\n" << endl;
volScalarField phiIB
(
  IOobject
  (
    "phiIB",
    runTime.timeName(),
    mesh,
    IOobject::READ_IF_PRESENT,
    IOobject::AUTO_WRITE
  ),
  mesh,
  dimensionedScalar("zero", dimensionSet(0, 2, -1, 0, 0), 0.0)
);

Info << "Reading field interFace\n" << endl;
volScalarField interFace
(
  IOobject
  (
    "interFace",
    runTime.timeName(),
    mesh,
    IOobject::READ_IF_PRESENT,
    IOobject::AUTO_WRITE
  ),
  mesh,
  dimensionedScalar("0", dimensionSet(0, 0, 0, 0, 0), 0.0)
);

Info << "Reading field refineMeshKeepStep\n" << endl;
volScalarField refineMeshKeepStep
(
  IOobject
  (
    "refineMeshKeepStep",
    runTime.timeName(),
    mesh,
    IOobject::READ_IF_PRESENT,
    IOobject::AUTO_WRITE
  ),
  mesh,
  dimensionedScalar("0", dimensionSet(0, 0, 0, 0, 0), 0.0)
);

Info << "Reading lagrangeMultiplierForce field\n" << endl;
volVectorField lmpf
(
  IOobject
  (
    "lmpf",
    runTime.timeName(),
    mesh,
    IOobject::NO_READ,
    IOobject::AUTO_WRITE
  ),
  mesh,
  dimensionedVector("zero", dimensionSet(0, 1, -2, 0, 0), vector(0, 0, 0))
);

// Info << "Reading field phiLmpf\n" << endl;
// volScalarField phiLmpf
// (
//   IOobject
//   (
//     "phiLmpf",
//     runTime.timeName(),
//     mesh,
//     IOobject::NO_READ,
//     IOobject::AUTO_WRITE
//   ),
//   mesh,
//   dimensionedScalar("zero", dimensionSet(0, 2, -2, 0, 0), 0.0)
// );
