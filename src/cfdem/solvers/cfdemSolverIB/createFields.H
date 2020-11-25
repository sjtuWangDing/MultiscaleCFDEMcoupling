Info << "Reading field p\n" << endl;
volScalarField p(
  IOobject(
    "p",
    runTime.timeName(),
    mesh,
    IOobject::MUST_READ,
    IOobject::AUTO_WRITE
  ),
  mesh
);

Info << "Reading physical velocity field U\n" << endl;
volVectorField U(
  IOobject(
    "U",
    runTime.timeName(),
    mesh,
    IOobject::MUST_READ,
    IOobject::AUTO_WRITE
  ),
  mesh
);

Info << "Reading volumeFraction field volumeFraction = (Vgas / Vparticle)\n" << endl;
volScalarField volumeFraction(
  IOobject(
    "volumeFraction",
    runTime.timeName(),
    mesh,
    IOobject::READ_IF_PRESENT,
    IOobject::AUTO_WRITE
  ),
  mesh,
  dimensionedScalar("zero", dimensionSet(0, 0, 0, 0, 0), 1.0)
);

Info << "Reading field interface\n" << endl;
volScalarField interface(
  IOobject(
    "interface",
    runTime.timeName(),
    mesh,
    IOobject::READ_IF_PRESENT,
    IOobject::AUTO_WRITE
  ),
  mesh,
  dimensionedScalar("0", dimensionSet(0, 0, 0, 0, 0), 0.0)
);

Info << "Reading field refineMeshKeepStep\n" << endl;
volScalarField refineMeshKeepStep(
  IOobject(
    "refineMeshKeepStep",
    runTime.timeName(),
    mesh,
    IOobject::READ_IF_PRESENT,
    IOobject::AUTO_WRITE
  ),
  mesh,
  dimensionedScalar("0", dimensionSet(0, 0, 0, 0, 0), 0.0)
);
