forAll(solidRegions, i)
{
    Info << "OF sending solid " << i << " receiving mesh to Elmer" << endl;
    solidReceiving.set
    (
        i,
        new Elmer<fvMesh>
        (
            solidRegions[i],
            -1
        )    
    );
}
