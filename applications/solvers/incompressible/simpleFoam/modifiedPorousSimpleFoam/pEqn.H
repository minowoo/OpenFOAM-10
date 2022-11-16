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
surfaceScalarField phiHbyA("phiHbyA", gammaf*fvc::flux(HbyA));

MRF.makeRelative(phiHbyA);

adjustPhi(phiHbyA, U, p);

while (simple.correctNonOrthogonal())
{
    tmp<fvScalarMatrix> tpEqn;

    if (pressureImplicitPorosity)
    {
        tmp<volTensorField> gammaTU;
        gammaTU = gamma*trTU();
        
        tpEqn = (fvm::laplacian(gammaTU(), p) == fvc::div(phiHbyA));
    }
    else
    {
        tmp<volScalarField> gammaAU;
        gammaAU = gamma*trAU();
            
        tpEqn = (fvm::laplacian(gammaAU(), p) == fvc::div(phiHbyA));
    }

    fvScalarMatrix pEqn = tpEqn.ref();

    pEqn.setReference
    (
        pressureReference.refCell(),
        pressureReference.refValue()
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
fvConstraints.constrain(U);
