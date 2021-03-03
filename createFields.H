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
// Reading heat of fusion
dimensionedScalar L(smProperties.lookup("L"));
// Reading solid heat capacity
dimensionedScalar cpS(smProperties.lookup("cp"));
// Reading thermal conductivity
dimensionedScalar lambdaS(smProperties.lookup("lambdaS"));
// Reading melting temperature
dimensionedScalar Tmelt(smProperties.lookup("Tmelt"));
// Reading mushy region momentum sink coefficient default is 1e5
dimensionedScalar Cu(smProperties.lookup("Cu"));
// Reading porosity coefficient
dimensionedScalar qD(smProperties.lookup("qD"));
// Reading relaxation coefficient 0.9
dimensionedScalar relaxAlpha(smProperties.lookup("relaxAlpha"));
// Reading reference density
dimensionedScalar rhoRef(smProperties.lookup("rhoRef"));
// Reading beta coeff
dimensionedScalar beta(smProperties.lookup("beta"));
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
