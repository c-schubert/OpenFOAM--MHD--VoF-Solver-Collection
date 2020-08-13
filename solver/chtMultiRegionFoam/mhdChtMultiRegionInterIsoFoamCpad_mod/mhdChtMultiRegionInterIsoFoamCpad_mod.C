/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2018 OpenCFD Ltd.
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Application
    chtMultiRegionTwoPhaseEulerFoam

Group
    grpHeatTransferSolvers

Description
    Transient solver for buoyant, turbulent fluid flow and solid heat
    conduction with conjugate heat transfer between solid and fluid regions.

    It solves a two-phase Euler approach on the fluid region.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"


#include "compressibleInterPhaseTransportModel.H"
#include "compressibleCourantNo.H"
#include "pimpleControl.H"
#include "fixedGradientFvPatchFields.H"
#include "regionProperties.H"
#include "solidRegionDiffNo.H"
#include "solidThermo.H"
#include "radiationModel.H"
#include "fvOptions.H"
#include "fvcSmooth.H"
#include "coordinateSystem.H"
#include "loopControl.H"

#include "advectionSchemes.H"
#include "twoPhaseMixtureThermo.H"

// needed?
#include "turbulentFluidThermoModel.H"
#include "pressureControl.H"

// Needed here?
#include "CMULES.H"
#include "EulerDdtScheme.H"
#include "localEulerDdtScheme.H"
#include "CrankNicolsonDdtScheme.H"
#include "subCycle.H"

#include "Elmer.H"

#include "cpad.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    // dynamicMesh features of compressivleInterFlow removed

    argList::addNote
    (
        "Transient solver for two compressible, non-isothermal immiscible fluids"
        "using VOF phase-fraction based interface capturing, with"
        "solid heat conduction with conjugate heat transfer "
        "between solid and fluid regions."
    );

    #define NO_CONTROL
    #define CREATE_MESH createMeshesPostProcess.H
    #include "postProcess.H"

    // #include "setRootCase.H"
    #include "setRootCaseLists.H"
    #include "createTime.H"

    #include "createMeshes.H"

    #include "createFields.H"
    #include "fluidCreateUf.H"

    #include "initContinuityErrs.H"
    #include "createTimeControls.H"
    #include "readSolidTimeControls.H"

    #include "compressibleMultiRegionCourantNo.H"
    #include "solidRegionDiffusionNo.H"
    #include "setInitialMultiRegionDeltaT.H"

    #include "createMhdFields.H"

    Info<< "\nStarting time loop\n" << endl;

    while (runTime.run())
    {
         bool updateElmerSolution = false;
         
        #include "readTimeControls.H"
        #include "readSolidTimeControls.H"
        #include "readPIMPLEControls.H"

        #include "compressibleMultiRegionCourantNo.H"
        #include "alphaMultiRegionCourantNo.H"
        #include "solidRegionDiffusionNo.H"
        #include "setMultiRegionDeltaT.H"

        ++runTime;

        Info<< "Time = " << runTime.timeName() << nl << endl;

        if (nOuterCorr != 1)
        {
            forAll(fluidRegions, i)
            {
                #include "storeOldFluidFields.H"
            }
        }

        // --- PIMPLE loop
        for (int oCorr=0; oCorr<nOuterCorr; ++oCorr)
        {
            const bool finalIter = (oCorr == nOuterCorr-1);
            const bool firstIter = (oCorr == 0);

            // needed or dynamic mesh relict?
            if (firstIter)
            {
                forAll(fluidRegions, i)
                {
                    #include "setRegionFluidFields.H" 
                    advector->surf().reconstruct();
                }
            }

            forAll(fluidRegions, i)
            {
                Info<< "\nSolving for fluid region "
                    << fluidRegions[i].name() << endl;
                #include "setRegionFluidFields.H"
                #include "setRegionMhdFluidFields.H"
                #include "readFluidMultiRegionPIMPLEControls.H" 
                #include "solveMhdFluid.H"
            }

            forAll(solidRegions, i)
            {
                Info<< "\nSolving for solid region "
                    << solidRegions[i].name() << endl;
                #include "setRegionSolidFields.H"
                #include "setRegionMhdSolidFields.H"
                #include "readSolidMultiRegionPIMPLEControls.H"
                #include "solveMhdSolid.H"
            }

            // Additional loops for energy solution only
            if (!oCorr && nOuterCorr > 1)
            {
                loopControl looping(runTime, pimple, "energyCoupling");

                while (looping.loop())
                {
                    Info<< nl << looping << nl;

                    forAll(fluidRegions, i)
                    {
                        Info<< "\nSolving for fluid region "
                            << fluidRegions[i].name() << endl;
                       #include "setRegionFluidFields.H"
                       #include "setRegionMhdFluidFields.H"
                       #include "readFluidMultiRegionPIMPLEControls.H"
                       frozenFlow = true;
                       #include "solveMhdFluid.H"
                    }

                    forAll(solidRegions, i)
                    {
                        Info<< "\nSolving for solid region "
                            << solidRegions[i].name() << endl;
                        #include "setRegionSolidFields.H"
                        #include "setRegionMhdSolidFields.H"
                        #include "readSolidMultiRegionPIMPLEControls.H"
                        #include "solveMhdSolid.H"
                    }
                }
            }
        } // End PIMPLE Loop

        // in case alpha1 - alpha1_old > maxRelDiff (in any fluid) -> update Elmer Solution
        forAll(fluidRegions, i)
        {
            #include "setRegionFluidFields.H"
            rho = alpha1*rho1 + alpha2*rho2;

            // Correct p_rgh for consistency with p and the updated densities
            p_rgh = p - rho*gh;
            p_rgh.correctBoundaryConditions();

            #include "setRegionMhdFluidFields.H"

            // #include "sendUpdateRegionMhdFluidFields.H"
            // modified content { 
            scalar maxRelDiff_local = (max(mag(alpha1_old - alpha1))).value();

            if
            ( 
                maxRelDiff_local>maxRelDiff 
                && (maxRelDiff<SMALL || maxRelDiff+SMALL<=1.0)
            ) 
            {
                updateElmerSolution = true;
                Info << "maxRelDifflocal = " << maxRelDiff_local << endl;
            }

            if(updateElmerSolution || !runTime.run()) 
            {
                Info << "FOAM: Continue Elmer!" << endl;
                alpha1_old = alpha1;


                alpha1f = min(max(alpha1, scalar(0)), scalar(1));
                elcond = alpha1f * elcond1 + (1 - alpha1f) * elcond2;
                // Send fields to Elmer

                forAll(alpha1,c)
                {
                    if (alpha1[c] < 0.1)
                    {
                        elcond[c] = elcond2.value();
                    }
                    if (elcond[c] < 0)
                    {
                        elcond[c] = 0;
                        Info << "elcond[c] smaller zero, this should not happen!" << endl;
                    }
                }

                sending.sendStatus(runTime.run());
                sending.sendScalar(elcond);

                // Receive fields from Elmer
                receiving.sendStatus(runTime.run());
                receiving.recvVector(JxB_recv);
                receiving.recvScalar(JH_recv);
            }

             // modified content }
        }

        forAll(solidRegions, i)
        {
            #include "setRegionSolidFields.H"
            #include "setRegionMhdSolidFields.H"
            #include "updateRegionMhdSolidFields.H"
        }

        runTime.write();

        forAll(fluidRegions, i)
        {
        // todo: modifications necessary if more that 1 fluid ...
        #include "setRegionFluidFields.H"
        #include "getLocalToGlobalFaceArray.h"

        FoamCpad::cpad cpas(mixture, runTime);
        cpas.appendTimeStepReport(
                                    runTime.value(),
                                    procFaceToGlobalFaceList
                                );
        }


        runTime.printExecutionTime(Info);
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
