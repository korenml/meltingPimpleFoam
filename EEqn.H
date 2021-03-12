{
for (int tCorr=0; tCorr<nCorrectorT; tCorr++)
{

    volScalarField& he = thermo.he();
	volScalarField alphaLNew = alphaL;
    volScalarField alphaL0 = alphaL;	
	//const volScalarField& T_ = thermo.T();
    fvScalarMatrix EEqn
    (
        fvm::ddt(rho, he) + fvm::div(phi, he)
      + fvc::ddt(rho, K) + fvc::div(phi, K)
      + (
            he.name() == "e"
          ? fvc::div
            (
                fvc::absolute(phi/fvc::interpolate(rho), U),
                p,
                "div(phiv,p)"
            )
          : -dpdt
        )
      - fvm::laplacian(turbulence->alphaEff(), he)
     ==
	  - L*fvc::ddt(rho, alphaL)
	  - L*fvc::div(phi, alphaL)
      +  rho*(U&g)
      + radiation->Sh(thermo, he)
      + fvOptions(rho, he)
    );

    EEqn.relax();

    fvOptions.constrain(EEqn);

    EEqn.solve();

    fvOptions.correct(he);

    thermo.correct();
    radiation->correct();

	alphaLNew = alphaL + relaxAlpha*cpS/L*(thermo.T()-Tmelt);
	alphaL = max(scalar(0), min(scalar(1), alphaLNew));
	alphaL.correctBoundaryConditions();

    // Check for phase fraction convergence
    scalarField residualField
    (
        mag(alphaL.internalField() - alphaL0.internalField())
    );

    scalar residual(gMax(residualField));

    Info<< "THERMO: iteration " << (tCorr + 1)
        << " residual: " << residual  << endl;
    if (residual < residualT)
    {
        break;
    }
}



}
