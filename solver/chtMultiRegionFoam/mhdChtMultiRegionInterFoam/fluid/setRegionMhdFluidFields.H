Elmer<fvMesh>& receiving  = fluidReceiving[i];
Elmer<fvMesh>& sending  = fluidSending[i];

volVectorField& JxB_recv = fluidJxB_recv[i];
volVectorField& JxB = fluidJxB[i];

volScalarField& JH_recv = fluidJH_recv[i];
volScalarField& JH = fluidJH[i];

dimensionedScalar elcond1 = fluidElCond1[i];
dimensionedScalar elcond2 = fluidElCond2[i];
volScalarField& elcond = fluidElCond[i];

volScalarField &alpha1f = fluidAlpha1f[i];
volScalarField& alpha1_old = fluidAlpha1_old[i];

scalar& maxRelDiff = fluidMaxRelDiff[i];


//JxB = JxB_recv*alpha1; // Korrektur (Levitating droplet ...)
JxB.field() = JxB_recv.field(); // New // per vol or absolute force?

JH.field() = JH_recv.field() ;

forAll(JH,c)
{
    if (JH[c]<0)
    {
        // Limit JH due to numerical errors in Elmerfem ...
        JH[c] = 0;
    }
}
