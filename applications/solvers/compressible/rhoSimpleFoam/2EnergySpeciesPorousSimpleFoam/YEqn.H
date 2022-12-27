{ 
    // Diffusivity of NO-air
    volScalarField Dmix = 1.43e-5 * Foam::pow(thermo.T(),1.75)/Foam::pow(293,1.75); 
    Dmix.dimensions().reset(dimensionSet(0, 2, -1, 0, 0, 0, 0));

    volScalarField km = Nu*Dmix/Dh;

    volScalarField ygSource = specificSurfaceArea*rho*km*(YG-YS);

    forAll(ygSource, celli)
    {
        if(gamma[celli]<1.0)
            {
                ygSource[celli] = ygSource[celli]; 
            } else
            {
                ygSource[celli] = ygSource[celli]*0.0; 
            }
    }		    

    fvScalarMatrix YGEqn
    (
       fvm::div(phi,YG)
    ==
       fvm::laplacian(gamma*rho*Dmix,YG)
       - ygSource
    );
    YGEqn.relax();
    YGEqn.solve();
    YG.max(0.0);
    //YG.min(1.0);

    volScalarField dimLessTs = ts;
    dimLessTs.dimensions().reset(dimensionSet(0,0,0,0,0,0,0));
    volScalarField reactionRate = 1.e4 * Foam::exp(-50000./8.314/dimLessTs);
    reactionRate = reactionRate*30/1000 * 0.1;  //unit mol/m2-s to kg/m2-s

//    volScalarField reactionRate = 4.6159e-8 * Foam::exp(1.3686e-2*dimLessTs)*0.5;
    reactionRate.dimensions().reset(dimensionSet(1,-2,-1,0,0,0,0));


    forAll(reactionRate, celli)
    {
        if(gamma[celli]<1.0)
    	{
	    reactionRate[celli] = reactionRate[celli];
    	} else
    	{
	    reactionRate[celli] = reactionRate[celli]*0.0;	
	}
    }		    
    //Info<<reactionRate<<endl;

    //YS = (-(1-gamma)*reactionRate + rho*km*YG)/(rho*km); 
    YS = (-reactionRate + rho*km*YG)/(rho*km); 
    fvConstraints.constrain(YS); 
    YS.max(0.0);
}
