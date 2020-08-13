
forAll(fluidRegions, i)
{
    Info << "OF sending  fluid " << i << " receiving mesh to Elmer" << endl;
    fluidReceiving.set
    (
        i,
        new Elmer<fvMesh>
        (
            fluidRegions[i],
            -1
        )    
    );
    
}
