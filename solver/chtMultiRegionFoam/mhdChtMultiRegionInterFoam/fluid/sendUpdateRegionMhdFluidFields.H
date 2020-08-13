scalar maxRelDiff_local = (max(mag(alpha1_old - alpha1))).value();

if
( 
    maxRelDiff_local>maxRelDiff 
    && (maxRelDiff<SMALL || maxRelDiff+SMALL<=1.0)
) 
{
    updateElmerSolution = true;
    Info << "maxRelDifflocal = " << maxRelDiff_local << endl;
}

if(updateElmerSolution || !runTime.run()) 
{
    Info << "FOAM: Continue Elmer!" << endl;
    alpha1_old = alpha1;

    alpha1f = min(max(alpha1, scalar(0)), scalar(1));
    elcond = alpha1f * elcond1 + (1 - alpha1f) * elcond2;
    // Send fields to Elmer

    sending.sendStatus(runTime.run());
    sending.sendScalar(elcond);

    // Receive fields from Elmer
    receiving.sendStatus(runTime.run());
    receiving.recvVector(JxB_recv);
    receiving.recvScalar(JH_recv);
}
