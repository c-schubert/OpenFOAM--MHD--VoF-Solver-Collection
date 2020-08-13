PtrList<Elmer<fvMesh>> fluidSending(fluidRegions.size());
PtrList<Elmer<fvMesh>> fluidReceiving(fluidRegions.size());

PtrList<volVectorField> fluidJxB_recv(fluidRegions.size());
PtrList<volVectorField> fluidJxB(fluidRegions.size());
PtrList<volScalarField> fluidJH_recv(fluidRegions.size());
PtrList<volScalarField> fluidJH(fluidRegions.size());

PtrList<volScalarField>  fluidAlpha1f(fluidRegions.size());

PtrList<scalar> fluidMaxRelDiff(fluidRegions.size());
PtrList<dimensionedScalar> fluidElCond1(fluidRegions.size());
PtrList<dimensionedScalar> fluidElCond2(fluidRegions.size());
PtrList<volScalarField> fluidElCond(fluidRegions.size());

PtrList<volScalarField> fluidAlpha1_old(fluidRegions.size());


forAll(fluidRegions, i)
{
     #include "setRegionFluidFields.H" // some things prob missing

    Info<< "Initializing field JxB_recv for fluid" << i << "\n" << endl;
    fluidJxB_recv.set
    (
        i,
        new volVectorField 
        (
            IOobject
            (
                "JxB_recv",
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh,
            dimensionedVector("zero", dimless,Foam::vector(0,0,0))
        )
    );

    Info<< "Initializing field JxB for fluid" << i << "\n" << endl;
    
    // Check!
    // Source term (UEqn) musst be in N per m³ - or N?
    // But is JxB sent by elmer allready N per m³ ? or N = kg m s^-2?
    // N per m³  =  kg m^-2 s^-2
    fluidJxB.set
    (
        i,
        new volVectorField 
        (
            IOobject
            (
                "JxB",
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            mesh,
            dimensionedVector("zero", dimensionSet(1, -2, -2, 0, 0, 0, 0),Foam::vector(0,0,0))
        )
    );

    
    Info<< "Initializing field JH_recv for fluid" << i << "\n" << endl;
    fluidJH_recv.set
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
    
    // W/m³ -> hoffentlich ...
    Info<< "Initializing field JH for fluid" << i << "\n" << endl;
    fluidJH.set
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

    fluidElCond1.set
    (
        i,
        new dimensionedScalar
        (
            "elcond", 
            dimensionSet(-1, -3, 3, 0, 0, 2, 0), 
            mixture.thermo1()
        )
    );

    Info << "fluid " << i << " mixture.thermo1 elcond: " << fluidElCond1[i] << endl;

    fluidElCond2.set
    (
        i,
        new dimensionedScalar
        (
            "elcond", 
            dimensionSet(-1, -3, 3, 0, 0, 2, 0), 
            mixture.thermo2()
        )
    );

    Info << "fluid " << i << " mixture.thermo2 elcond: " << fluidElCond2[i] << endl;

    // Maximum allowable relative field difference in any cell comparing
    // to last Elmer update. This option controls frequency Elmer is called.
    //     1     - once
    //     (0,1) - when relative difference in any cell exceeds given value
    //     0     - every iteration

    fluidMaxRelDiff.set
    (
        i,
        new scalar
        (
            readScalar(mixture.lookup("maxRelDiff"))
        )
    );

    
    
    Info<< "Initializing field fluidAlpha1_old\n" << endl;
    fluidAlpha1_old.set
    (
        i,
        new volScalarField
        (
            IOobject
            (
                "alpha_old",
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh,
            dimensionedScalar("zero", dimless, 0)
        )
    );

    Info<< "Initializing field alpha1f (safty contains non negative alpha values ...)\n" << endl;
    fluidAlpha1f.set
    (
        i,
        new volScalarField
        (
            IOobject
            (
                "alpha1f",
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            alpha1
        )
    );


    Info<< "Initializing fluid " << i << " field ElCond\n" << endl;
    fluidElCond.set
    (
        i,
        new volScalarField
        (
            IOobject
            (
                "elcond",
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE // NO_WRITE
            ),
            mesh,
            dimensionedScalar("zero", dimensionSet(-1, -3, 3, 0, 0, 2, 0), 0)
        )
    );

    Info << "Max alpha1 init " << gMax(alpha1) << endl;
    Info << "Min alpha1 init " << gMin(alpha1) << endl;

    volScalarField& alpha1f = fluidAlpha1f[i];
    alpha1f = min(max(alpha1, scalar(0)), scalar(1));

    Info << "Max alpha1f init " << gMax(alpha1f) << endl;
    Info << "Min alpha1f init " << gMin(alpha1f) << endl;

    fluidElCond[i] = alpha1f *  fluidElCond1[i] + (1 - alpha1f) * fluidElCond2[i];
    fluidElCond[i] = max(fluidElCond[i], dimensionedScalar("zero",dimensionSet(-1, -3, 3, 0, 0, 2, 0), 0));

    Info << "Max elcond init " << gMax(fluidElCond[i]) << endl;
    Info << "Min elcond init " << gMin(fluidElCond[i]) << endl;

    Info << "Init fluid " << i << " finished!" << endl;
}
