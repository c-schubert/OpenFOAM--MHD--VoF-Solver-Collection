forAll(fluidRegions, i)
{
    fluidSending[i].sendStatus(1);
    Info << "OF sending fluid " << i << " elcond to Elmer" << endl;
    fluidSending[i].sendScalar(fluidElCond[i]);
}
