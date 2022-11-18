{ 
    // for alkane hydrocarbons
    volScalarField Dmix = 2.1e-5 * Foam::pow(thermo.T(),1.75)/Foam::pow(293,1.75); //methane // m2/s
    Dmix.dimensions().reset(dimensionSet(0, 2, -1, 0, 0, 0, 0));

    scalar carbon_number = 40; // tetracontane
    Dmix = -0.02198 * Dmix * carbon_number + 1.021978 * Dmix;

    fvScalarMatrix YGEqn
    (
       fvm::div(phi,YG)
    ==
       fvm::laplacian(gamma*rho*Dmix,YG)
    );
    YGEqn.relax();
    YGEqn.solve();
    YG.max(0.0);
    //YG.min(1.0);

    /*
    volScalarField TdegC = thermo.T();
    TdegC.dimensions().reset(dimensionSet(0,0,0,0,0,0,0));

    volScalarField psat = exp (7.16731 - 2422.3/(TdegC-273+40.1)); //mmHg 
    psat.dimensions().reset(dimensionSet(1,-1,-2,0,0,0,0));
    S = p/(psat*133.33);
    */
}
