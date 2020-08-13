// not needed (maybe in the future for T-prop coupling ...)
// PtrList<Elmer<fvMesh>> solidSending(solidRegions.size());
PtrList<Elmer<fvMesh>> solidReceiving(solidRegions.size());
PtrList<volScalarField> solidJH_recv(solidRegions.size());
PtrList<volScalarField> solidJH(solidRegions.size());

forAll(solidRegions, i)
{
     #include "setRegionSolidFields.H" // some things prob missing

    Info<< "Initializing field JH_recv for solid " << i << "\n" << endl;
    solidJH_recv.set
    (
        i,
        new volScalarField
        (
            IOobject
            (
            "JH_recv",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
            ),
            mesh,
            dimensionedScalar("zero", dimless, 0)
        )
    );


    Info<< "Initializing field JH for solid " << i << "\n" << endl;
    solidJH.set
    (
        i,
        new volScalarField
        (
            IOobject
            (
                "JH",
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            mesh,
            dimensionedScalar("zero", dimensionSet(1, -1, -3, 0, 0, 0, 0), 0)
        )
    );

    Info << "Init solid " << i << " finished!" << endl;
}
