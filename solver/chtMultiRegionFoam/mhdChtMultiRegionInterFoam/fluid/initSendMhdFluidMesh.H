forAll(fluidRegions, i)
{
    Info << "OF sending fluid " << i << " sending mesh to Elmer" << endl;
    fluidSending.set
    (
        i,
        new Elmer<fvMesh>
        (
            fluidRegions[i],
            1
        )    
    );
}
