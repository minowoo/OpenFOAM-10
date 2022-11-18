const volScalarField& psi = thermo.psi();

tmp<volVectorField> tHbyA;
surfaceScalarField gammaf(fvc::interpolate(gamma));
if (pressureImplicitPorosity)
{
    tHbyA = constrainHbyA(trTU()&UEqn.H(), U, p);
}
else
{
    tHbyA = constrainHbyA(trAU()*UEqn.H(), U, p);
}
volVectorField& HbyA = tHbyA.ref();

tUEqn.clear();

bool closedVolume = false;

//surfaceScalarField phiHbyA("phiHbyA", fvc::flux(HbyA));
surfaceScalarField phiHbyA("phiHbyA", gammaf*fvc::flux(rho*HbyA));

MRF.makeRelative(fvc::interpolate(rho), phiHbyA);

closedVolume = adjustPhi(phiHbyA, U, p);

while (simple.correctNonOrthogonal())
{
    tmp<fvScalarMatrix> tpEqn;

    if (pressureImplicitPorosity)
    {
        tmp<volTensorField> gammaTU;
        gammaTU = gamma*rho*trTU();

        //tpEqn = (fvm::laplacian(trTU(), p) == fvc::div(phiHbyA));
        tpEqn = (fvm::laplacian(gammaTU(), p) 
                +fvOptions(psi, p, rho.name())
                == fvc::div(phiHbyA));
    }
    else
    {
        tmp<volScalarField> gammaAU;
        gammaAU = gamma*rho*trAU();

        //tpEqn = (fvm::laplacian(trAU(), p) == fvc::div(phiHbyA));
        tpEqn = (fvm::laplacian(gammaAU(), p) 
                 + fvOptions(psi, p, rho.name())
                == fvc::div(phiHbyA));
    }

    fvScalarMatrix pEqn = tpEqn.ref();

        pEqn.setReference
        (
            pressureControl.refCell(),
            pressureControl.refValue()
        );

    pEqn.solve();

    if (simple.finalNonOrthogonalIter())
    {
        phi = phiHbyA - pEqn.flux();
    }
}

#include "continuityErrs.H"

// Explicitly relax pressure for momentum corrector
p.relax();

if (pressureImplicitPorosity)
{
    U = HbyA - (trTU()&fvc::grad(p));
}
else
{
    U = HbyA - (trAU()*fvc::grad(p));
}

U.correctBoundaryConditions();
fvOptions.correct(U);

    pressureControl.limit(p);

    // For closed-volume cases adjust the pressure and density levels
    // to obey overall mass continuity
    if (closedVolume)
    {
        p += (initialMass - fvc::domainIntegrate(psi*p))
            /fvc::domainIntegrate(psi);
    }


rho = thermo.rho();
rho = max(rho, rhoMin);
rho = min(rho, rhoMax);
rho.relax();
Info<< "rho max.min : " << max(rho).value() << " " << min(rho).value() << endl;
