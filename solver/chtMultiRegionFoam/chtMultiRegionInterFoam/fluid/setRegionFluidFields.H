    fvMesh& mesh = fluidRegions[i];

    twoPhaseMixtureThermo& mixture  = mixtureFluid[i];

    volScalarField& alpha1(mixture.alpha1());
    volScalarField& alpha2(mixture.alpha2());

    volScalarField& rho1 = mixture.thermo1().rho();
    volScalarField& rho2 = mixture.thermo2().rho();

    volVectorField& U = UFluid[i];
    surfaceScalarField& phi = phiFluid[i];
    volScalarField& rho = rhoFluid[i];

    IOMRFZoneList& MRF = MRFfluid[i];
    fv::options& fvOptions = fvOptionsFluid[i];

    volScalarField& p_rgh = p_rghFluid[i];

    surfaceScalarField& rhoPhi = rhoPhiFluid[i];

    volScalarField& p = mixture.p();
    volScalarField& T = mixture.T();
    const volScalarField& psi1 = mixture.thermo1().psi();
    const volScalarField& psi2 = mixture.thermo2().psi();

    //const dimensionedVector& g = gFluid[i];
    const volScalarField& gh = ghFluid[i];
    const surfaceScalarField& ghf = ghfFluid[i];

    const bool alphaRestart = alphaRestartFluid[i];
    volScalarField& dgdt = dgdtFluid[i];
    surfaceScalarField& alphaPhi10 = alphaPhi10Fluid[i];

    // MULES Correction
    tmp<surfaceScalarField> talphaPhi1Corr0;

    compressibleInterPhaseTransportModel&  turbulence = turbulenceFluid[i];
    
    volScalarField& K = KFluid[i];

    const dimensionedScalar initialMass
    (
        "initialMass",
        dimMass,
        initialMassFluid[i]
    );

    bool frozenFlow = frozenFlowFluid[i];

    //const label pRefCell = pRefCellFluid[i];
    //const scalar pRefValue = pRefValueFluid[i];

    // pimpleControl& pimple = pimpleFluid[i];

    const dimensionedScalar& pMin = pMinFluid[i];

    // Debug fields
/*
    volScalarField& faceRegimes = faceRegimesFluid[i];
    volScalarField& qc = qcFluid[i];
    volScalarField& qFilm = qFilmFluid[i];
    volScalarField& htcFilmBoiling = htcFilmBoilingFluid[i];
    volScalarField& qtb = qtbFluid[i];
    volScalarField& qSubCool = qSubCoolFluid[i];
    volScalarField& CHFtotal = CHFtotalFluid[i];
    volScalarField& Tdnb = TdnbFluid[i];
    volScalarField& phiTb = phiFluid[i];
*/
