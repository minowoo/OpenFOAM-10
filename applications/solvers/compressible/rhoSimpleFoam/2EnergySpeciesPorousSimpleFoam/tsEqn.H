{
    const dimensionedScalar coeffA
    ("coeffA", dimensionSet(0, 0, 0, 1, 0, 0, 0), 0.9558);

    volScalarField kappaSc = coeffA - 2.09e-4 * ts; //ks for ceramic

    volScalarField dimLessTs = ts;
    dimLessTs.dimensions().reset(dimensionSet(0,0,0,0,0,0,0));

    volScalarField kappaSm = -1.4087 
	                    +1.3982 * Foam::log10(dimLessTs) // ks for metal 
	                    +0.2543 * Foam::pow(Foam::log10(dimLessTs),2) 
	                    -0.626  * Foam::pow(Foam::log10(dimLessTs),3) 
	                    +0.2334 * Foam::pow(Foam::log10(dimLessTs),4) 
	                    +0.4256 * Foam::pow(Foam::log10(dimLessTs),5) 
	                    -0.4658 * Foam::pow(Foam::log10(dimLessTs),6) 
	                    +0.165  * Foam::pow(Foam::log10(dimLessTs),7) 
	                    -0.0199 * Foam::pow(Foam::log10(dimLessTs),8); 
    kappaSm = Foam::pow(10,kappaSm);

    kappaSm.dimensions().reset(dimensionSet(0,0,0,1,0,0,0));

    volScalarField kappaS = kappaSm;

    forAll(kappaS,celli)
    {
	if(gamma[celli]<0.7)
	{
	    kappaS[celli] = kappaSc[celli];
	} 
	    //Info<<kappaSc[celli]<<","<<kappaSm[celli]<<endl;
    }

    kappaS.dimensions().reset(dimensionSet(1,1,-3,-1,0,0,0));

    volScalarField kappa = thermo.kappa();
    volScalarField h_int = Nu*kappa/Dh;

    volScalarField tsSource = specificSurfaceArea*h_int*(ts-thermo.T());

    forAll(tsSource, celli)
    {
	if(gamma[celli]<1.0)
	{
            tsSource[celli] = tsSource[celli]; // /(1.0-gamma[celli]); //*mesh.V()[celli];
	} else
	{
	    tsSource[celli] = tsSource[celli]*0.0; //mesh.V()[celli]*0.0;
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
