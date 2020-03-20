Info << "Reading transportProperties\n" << endl;
IOdictionary transportProperties
(
  IOobject
  (
    "transportProperties",
    runTime.constant(),
    mesh,
    IOobject::MUST_READ_IF_MODIFIED,
    IOobject::NO_WRITE
  )
);

dimensionedScalar nu
(
  "nu",
  dimViscosity,
  transportProperties.lookup("nu")
);

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


Info << "Reading field U\n" << endl;
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

// #include "createPhi.H"
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
  fvc::flux(U)
);

label pRefCell = 0;
scalar pRefValue = 0.0;
setRefCell(p, mesh.solutionDict().subDict("PISO"), pRefCell, pRefValue);
mesh.setFluxRequired(p.name());
