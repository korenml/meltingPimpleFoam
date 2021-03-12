Info<< "Reading thermophysical properties\n" << endl;

autoPtr<rhoThermo> pThermo(rhoThermo::New(mesh));
rhoThermo& thermo = pThermo();
thermo.validate(args.executable(), "h", "e");

volScalarField rho
(
    IOobject
    (
        "rho",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    thermo.rho()
);

volScalarField h
(
    IOobject
    (
        "h",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    thermo.he()
);

volScalarField& p = thermo.p();

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

#include "compressibleCreatePhi.H"

pressureControl pressureControl(p, rho, pimple.dict(), false);

Info<< "Creating turbulence model\n" << endl;
autoPtr<compressible::turbulenceModel> turbulence
(
    compressible::turbulenceModel::New
    (
        rho,
        U,
        phi,
        thermo
    )
);


#include "readGravitationalAcceleration.H"
#include "readhRef.H"
#include "gh.H"


Info<< "Reading field p_rgh\n" << endl;
volScalarField p_rgh
(
    IOobject
    (
        "p_rgh",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

// Force p_rgh to be consistent with p
p_rgh = p - rho*gh;

mesh.setFluxRequired(p_rgh.name());

label pRefCell = 0;
scalar pRefValue = 0.0;

if (p_rgh.needReference())
{
    setRefCell
    (
        p,
        p_rgh,
        pimple.dict(),
        pRefCell,
        pRefValue
    );

    p += dimensionedScalar
    (
        "p",
        p.dimensions(),
        pRefValue - getRefCellValue(p, pRefCell)
    );
}

dimensionedScalar initialMass("initialMass", fvc::domainIntegrate(rho));

#include "createDpdt.H"

#include "createK.H"

#include "createMRF.H"
#include "createRadiationModel.H"
#include "createFvOptions.H"


const dimensionedScalar rhoMax("rhoMax", dimDensity, GREAT, pimple.dict());
const dimensionedScalar rhoMin("rhoMin", dimDensity, Zero, pimple.dict());

Info << "Reading liquid content" << endl;
volScalarField alphaL
(
    IOobject
    (
        "alphaL",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

Info << "Reading solidification melting properties" << endl;

IOdictionary smProperties
(
	IOobject
	(
		"smProperties",
		runTime.constant(),
		mesh,
		IOobject::MUST_READ,
		IOobject::NO_WRITE
	)
);

scalar nCorrectorT(smProperties.getOrDefault<scalar>("nCorrectorT", 100));
//int nCorrectorT = 10;
//label resT(smProperties.lookup("residualT"));
scalar residualT(smProperties.getOrDefault<scalar>("residualT", 1e-4));

// Reading heat of fusion
dimensionedScalar L("L", dimensionSet(0, 2, -2, 0, 0, 0, 0), smProperties.get<scalar>("L"));
// Reading solid heat capacity
dimensionedScalar cpS("cpS", dimensionSet(0, 2, -2, -1, 0, 0, 0), smProperties.get<scalar>("cpS"));
// Reading thermal conductivity
dimensionedScalar lambdaS("lambdaS", dimensionSet(1, 1, -3, -1, 0, 0, 0), smProperties.get<scalar>("lambdaS"));
// Reading melting temperature
dimensionedScalar Tmelt("Tmelt", dimensionSet(0, 0, 0, 1, 0, 0, 0), smProperties.get<scalar>("Tmelt"));
// Reading melting temperature difference
dimensionedScalar deltaT("deltaT", dimensionSet(0, 0, 0, 1, 0, 0, 0), smProperties.get<scalar>("deltaT"));
dimensionedScalar Tdim("Tdim", dimensionSet(0, 0, 0, 1, 0, 0, 0), smProperties.get<scalar>("Tdim"));
dimensionedScalar TS = Tmelt - 0.5*deltaT;
dimensionedScalar TL = Tmelt + 0.5*deltaT;
// Reading mushy region momentum sink coefficient default is 1e5
dimensionedScalar Cu("Cu", dimensionSet(0, 0, -1, 0, 0, 0, 0), smProperties.getOrDefault<scalar>("Cu", 1e5));
// Reading porosity coefficient
dimensionedScalar qD("qD", dimensionSet(0, 0, 0, 0, 0, 0, 0), smProperties.getOrDefault<scalar>("qD", 1e-3));
// Reading relaxation coefficient 0.9
dimensionedScalar relaxAlpha("relaxAlpha", dimensionSet(0, 0, 0, 0, 0, 0, 0), smProperties.getOrDefault<scalar>("relaxAlpha", 0.9));
// Creating momentum sink term
volScalarField DCu
(
    IOobject
    (
        "DCu",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    rho*Cu*pow(1.0-alphaL, 2)/(pow(alphaL, 3)+qD)
);
