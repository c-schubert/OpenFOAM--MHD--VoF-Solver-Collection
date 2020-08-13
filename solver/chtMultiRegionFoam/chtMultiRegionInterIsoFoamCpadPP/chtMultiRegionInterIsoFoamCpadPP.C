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
#include "fixedGradientFvPatchFields.H"
#include "regionProperties.H"
#include "fvOptions.H"
#include "fvcSmooth.H"
#include "coordinateSystem.H"

#include "advectionSchemes.H"
#include "twoPhaseMixtureThermo.H"
#include "turbulentFluidThermoModel.H"
#include "pressureControl.H"
#include "dynamicFvMesh.H"


#include "cpad.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    argList::addNote
    (
        "Cpad post process utility for chtMultiRegionInterFoam cases"
    );

    #define NO_CONTROL
    // #define CREATE_MESH createMeshesPostProcess.H
    // #include "postProcess.H"

    // #include "setRootCase.H"
    #include "setRootCaseLists.H"
    #include "createTime.H"

    regionProperties rp(runTime);
    #include "createFluidMeshes.H"
    #include "createFluidFields.H"
    #include "fluidCreateUf.H"


    Info<< "-----------------------------------------------------------------"
        << endl;
    Info<< "Starting continous phase area detection for processed cases:" 
        << endl;
    Info<< "\tCase directory: " << runTime.path() << endl;
    Info<< "\tCase name: " << runTime.caseName() << endl;
    Info<< "-----------------------------------------------------------------\n"
        << endl;

    Foam::instantList availTime(runTime.times());

    {
        forAll(availTime, t)
        {
            runTime.setTime(availTime[t].value(), 1);

            forAll(fluidRegions, i)
            {
                // todo: modifications necessary if more that 1 fluid ...
                #include "createFluidMeshes.H"
                #include "createFluidFields.H"
                #include "setRegionFluidFields.H"
                #include "getLocalToGlobalFaceArray.h"

                FoamCpad::cpad cpas(mixture, runTime);
                cpas.appendTimeStepReport
                (
                                runTime.value(),
                                procFaceToGlobalFaceList         
                );
            }


            runTime.printExecutionTime(Info);
        }
    }

    Info<< "Cpa detection finished\n" << endl;

    return 0;
}


// ************************************************************************* //
