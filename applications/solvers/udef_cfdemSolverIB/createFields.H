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
    IOobject::MUST_READ,
    IOobject::AUTO_WRITE
  ),
  mesh
);

Info << "\nCreating dummy density field rho = 1\n" << endl;
volScalarField rho
(
  IOobject
  (
    "rho",
    runTime.timeName(),
    mesh,
    IOobject::READ_IF_PRESENT,
    IOobject::AUTO_WRITE
  ),
  mesh,
  dimensionedScalar("0", dimensionSet(1, -3, 0, 0, 0), 1.0)
);

// 引入 phiIB 通量的原因：因为在速度映射过程中，在混合流体中（颗粒内部和外部），速度是连续的，但是在颗粒边界处，会出现颗粒不连续的情况，所以这里引入 phiIB 这个标量。
Info << "Reading field phiIB\n" << endl;
volScalarField phiIB
(
  IOobject
  (
    "phiIB",
    runTime.timeName(),
    mesh,
    IOobject::MUST_READ,
    IOobject::AUTO_WRITE
  ),
  mesh
);

Info << "\nReading voidfraction field voidfraction = (Vgas/Vparticle)\n" << endl;
volScalarField voidfraction
(
  IOobject
  (
    "voidfraction",
    runTime.timeName(),
    mesh,
    IOobject::MUST_READ,
    IOobject::AUTO_WRITE
  ),
  mesh
);

#include "createPhi.H"
// Info << "Reading/calculating face flux field phi\n" << endl;
// surfaceScalarField phi
// (
//   IOobject
//   (
//     "phi",
//     runTime.timeName(),
//     mesh,
//     IOobject::READ_IF_PRESENT,
//     IOobject::AUTO_WRITE
//   ),
//   fvc::flux(U)
// );

label pRefCell = 0;
scalar pRefValue = 0.0;
setRefCell(p, mesh.solutionDict().subDict("PISO"), pRefCell, pRefValue);

singlePhaseTransportModel laminarTransport(U, phi);

autoPtr<incompressible::turbulenceModel> turbulence
(
  incompressible::turbulenceModel::New(U, phi, laminarTransport)
);

// dynamic grid
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
