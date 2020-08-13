/*--------------------------------*- C++ -*----------------------------------*\
| =========                 |                                                 |
| \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox           |
|  \\    /   O peration     | Version:  v1912                                 |
|   \\  /    A nd           | Website:  www.openfoam.com                      |
|    \\/     M anipulation  |                                                 |
\*---------------------------------------------------------------------------*/
FoamFile
{
    version     2.0;
    format      ascii;
    class       dictionary;
    location    "constant/mixture";
    object      thermophysicalProperties;
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

thermoType
{
    type            heRhoThermo;
    mixture         pureMixture;
    transport       const;
    thermo          hConst;
    equationOfState rPolynomial;
    specie          specie;
    energy          sensibleInternalEnergy;
}

mixture
{
    specie
    {
        molWeight   55.691; // 1.2344
    }
    thermodynamics
    {
        Cp          820;
        Hf          0;
    }
    equationOfState
    {
        C    (1.3836E-4 -2.0193E-8 1.2537E-11 0 0);
    }
    transport
    {
        mu          0.0065;
        Pr          0.167; // kappa (lambda = 32)
    }
}

elcond [-1 -3  3 0 0 2 0] 714285.71;

// ************************************************************************* //
