    // from compressibleInterFoam (unchanged)
    MRF.correctBoundaryVelocity(U); // needed?
    
    fvVectorMatrix UEqn
    (
        fvm::ddt(rho, U) + fvm::div(rhoPhi, U)
      - fvm::Sp(contErr, U)
      + MRF.DDt(rho, U)
      + turbulence.divDevRhoReff(U)
     ==
        fvOptions(rho, U)
        + JxB
    );

    UEqn.relax();

    fvOptions.constrain(UEqn);

    if (momentumPredictor)
    {
        solve
        (
            UEqn
         ==
            fvc::reconstruct
            (
                (
                     mixture.surfaceTensionForce()
                    - ghf*fvc::snGrad(rho)
                    - fvc::snGrad(p_rgh)
                ) * mesh.magSf()
            )
        );

        fvOptions.correct(U);

        K = 0.5*magSqr(U);
    }
