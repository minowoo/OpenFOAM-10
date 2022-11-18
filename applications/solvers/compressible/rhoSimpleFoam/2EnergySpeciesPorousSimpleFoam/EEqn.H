{
    volScalarField& he = thermo.he();
    volScalarField kappa = thermo.kappa();
    volScalarField h_int = Nu*kappa/Dh;
    volScalarField hSource = specificSurfaceArea/gamma*h_int*(thermo.T()-ts);

    forAll(hSource, celli)
    {
	if(gamma[celli]<1.0)
	{	
	    hSource[celli] = hSource[celli]*mesh.V()[celli];
	} else
	{
	    hSource[celli] = hSource[celli]*mesh.V()[celli]*0.0;
	}
    }

    fvScalarMatrix EEqn
    (
        fvm::div(phi, he)
      + (
            he.name() == "e"
          ? fvc::div(phi, volScalarField("Ekp", 0.5*magSqr(U) + p/rho))
          : fvc::div(phi, volScalarField("K", 0.5*magSqr(U)))
        )
      + thermophysicalTransport->divq(he)
      + hSource
     ==
        fvModels.source(rho, he)
    );

    EEqn.relax();

    fvConstraints.constrain(EEqn);

    EEqn.solve();

    fvConstraints.constrain(he);

    thermo.correct();
}
