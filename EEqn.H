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
	  + L*(fvc::ddt(rho, alphaL))
     ==
        rho*(U&g)
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
	volScalarField& alphaLold = alphaL;
	alphaL = min(0, max(1, alphaL + relax*CpS*(T-Tmelt)/L));
	alphaL.correctBoundaryConditions();

	DCu = Cu*pow(1.0-alphaL, 2)/(pow(alphaL, 3)+q);
}
