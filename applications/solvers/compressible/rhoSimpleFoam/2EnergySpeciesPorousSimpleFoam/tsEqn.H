{
    const dimensionedScalar coeffA
    ("coeffA", dimensionSet(0, 0, 0, 1, 0, 0, 0), 0.9558);

    volScalarField kappaS = coeffA - 2.09e-4 * ts;
    kappaS.dimensions().reset(dimensionSet(1,1,-3,-1,0,0,0));

    volScalarField kappa = thermo.kappa();
    volScalarField h_int = Nu*kappa/Dh;
    volScalarField tsSource = specificSurfaceArea/gamma*h_int*(ts-thermo.T());

    forAll(tsSource, celli)
    {
	if(gamma[celli]<1.0)
	{
            tsSource[celli] = tsSource[celli]*mesh.V()[celli];
	} else
	{
	    tsSource[celli] = tsSource[celli]*mesh.V()[celli]*0.0;
	}
    }

    fvScalarMatrix tsEqn
    (
    	-fvm::laplacian(kappaS,ts)
	+ tsSource 
    );

    tsEqn.relax();

    tsEqn.solve();
}
