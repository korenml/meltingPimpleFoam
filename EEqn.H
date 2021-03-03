{
    volScalarField& he = thermo.he();
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
        rho*(U&g)
	  - L*(fvc::ddt(rho, alphaL))
      + radiation->Sh(thermo, he)
      + fvOptions(rho, he)
    );

    EEqn.relax();

    fvOptions.constrain(EEqn);

    EEqn.solve();

    fvOptions.correct(he);

    thermo.correct();
    radiation->correct();
	// Recalculate liquid content alpha
//	const volScalarField& T_ = thermo.T(); //mesh.lookupObject<volScalarField>("T");
	volScalarField alphaLnew =  alphaL + relaxAlpha*thermo.Cp()*(thermo.T()-Tmelt)/L;
	alphaL = max(scalar(0), min(scalar(1), alphaLnew));
	//alphaL.correctBoundaryConditions();

	//DCu = rho*Cu*pow(1.0-alphaL, 2)/(pow(alphaL, 3)+q);
}
