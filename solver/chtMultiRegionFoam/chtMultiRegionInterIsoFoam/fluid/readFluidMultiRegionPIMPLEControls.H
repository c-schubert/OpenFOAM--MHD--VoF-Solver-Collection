
const dictionary& pimpleDict = mesh.solutionDict().subDict("PIMPLE");

const int nCorr =
    pimpleDict.lookupOrDefault<int>("nCorrectors", 1);

const int nNonOrthCorr =
    pimpleDict.lookupOrDefault<int>("nNonOrthogonalCorrectors", 0);

const bool momentumPredictor =
    pimpleDict.lookupOrDefault("momentumPredictor", true);

// needed?
int nEnergyCorrectors = 
    pimpleDict.lookupOrDefault<int>("nEnergyCorrectors", 1);

const bool turbOnFinalIterOnly = 
    pimpleDict.lookupOrDefault<Switch>("turbOnFinalIterOnly", true);

const bool transonic = 
    mesh.solutionDict().lookupOrDefault<Switch>("transonic", false);
