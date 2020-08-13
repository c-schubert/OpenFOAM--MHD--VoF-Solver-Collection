    scalar CoNum = -GREAT;

    forAll(fluidRegions, regioni)
    {
        CoNum = max
        (
            Foam::compressibleCourantNo
            (
                fluidRegions[regioni],
                runTime,
                rhoFluid[regioni],
                phiFluid[regioni]
            ),
            CoNum
        );
    }
/*
    forAll(porousFluidRegions, porousi)
    {
        CoNum = max
        (
            compressibleCourantNo
            (
                porousFluidRegions[porousi],
                runTime,
                rhoPorous[porousi],
                phiPorous[porousi]
            ),
            CoNum
        );
    }
*/
