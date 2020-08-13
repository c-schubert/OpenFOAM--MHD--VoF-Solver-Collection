// from compressibleInterFoam (unchanged)


if (finalIter)
{
    mesh.data::add("finalIteration", true);
}


if (frozenFlow)
{   
    // from compressibleAlphaEqnSubCycle
    #include "contErr.H"
    volScalarField divUp("divUp", fvc::div(fvc::absolute(phi, U), p));
    #include "mhdTEqn.H"
}
else
{
    #include "alphaControls.H" // which version should be used -> gerneal icAlpha!
    #include "compressibleAlphaEqnSubCycle.H"

    turbulence.correctPhasePhi();

    #include "mhdUEqn.H"
    volScalarField divUp("divUp", fvc::div(fvc::absolute(phi, U), p));
    #include "mhdTEqn.H" // EEqn would be needed to include radiation
                    // T sufficient if no rad?

    // --- Pressure corrector loop
    // --- PISO loop
    for (int corr=0; corr<nCorr; corr++)
    {
        #include "pEqn.H"
    }

    if((turbOnFinalIterOnly && finalIter) || !turbOnFinalIterOnly)
    {
        turbulence.correct();
    }
}

if (finalIter)
{
    mesh.data::remove("finalIteration");
}
