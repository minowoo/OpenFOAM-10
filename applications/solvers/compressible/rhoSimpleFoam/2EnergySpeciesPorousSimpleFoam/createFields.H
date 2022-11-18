    Info<< "Reading thermophysical properties\n" << endl;

    autoPtr<fluidThermo> pThermo
    (
        fluidThermo::New(mesh)
    );
    fluidThermo& thermo = pThermo();
thermo.validate(args.executable(), "h", "e");

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
        thermo.rho()
    );

    volScalarField& p = thermo.p();

    Info<< "Reading field ts\n" << endl;
        volScalarField ts
        (
            IOobject
            (
                "ts",
                runTime.timeName(),
                mesh,
                IOobject::MUST_READ,
                IOobject::AUTO_WRITE
            ),
            mesh
        );


Info<< "Reading field U\n" << endl;
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

Info<< "Reading field gamma\n" << endl;
volScalarField gamma 
(
    IOobject
    (
        "gamma",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

Info<< "Reading field YG\n" << endl;
volScalarField YG
(
    IOobject
    (
        "YG",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

Info<< "Reading field S\n" << endl;
volScalarField S
(
    IOobject
    (
        "S",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

// #include "createPhi.H"

Info<< "Reading/calculating face flux field phi\n" << endl;
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
    fvc::interpolate(gamma)*linearInterpolate(rho*U) & mesh.Sf() //gammaPhi 
);

pressureReference pressureReference
(
     p,
     simple.dict(),
     thermo.incompressible()
);

mesh.schemes().setFluxRequired(p.name());

Info<< "Creating turbulence model\n" << endl;
autoPtr<compressible::momentumTransportModel> turbulence
(
     compressible::momentumTransportModel::New
     (
          rho,
          U,
          phi,
          thermo
     )
);

Info<< "Creating thermophysical transport model\n" << endl;
autoPtr<fluidThermophysicalTransportModel> thermophysicalTransport
(
     fluidThermophysicalTransportModel::New(turbulence(), thermo)
);

dimensionedScalar initialMass = fvc::domainIntegrate(rho);

#include "createMRF.H"
#include "createFvModels.H"
#include "createFvConstraints.H"

const dimensionedScalar Nu ("Nu", dimensionSet(0,0,0,0,0,0,0), 5); // Nusselt number
const dimensionedScalar Dh ("Dh", dimensionSet(0,1,0,0,0,0,0), 2e-3); // Hydraulic diameter of channel, m
const dimensionedScalar specificSurfaceArea
("specificSurfaceArea", dimensionSet(0,-1,0,0,0,0,0), 1000); // specific surfac area, m-1
