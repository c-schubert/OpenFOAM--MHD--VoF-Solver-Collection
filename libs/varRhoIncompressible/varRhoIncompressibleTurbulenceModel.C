/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2013-2014 OpenFOAM Foundation
     \\/     M anipulation  |
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

\*---------------------------------------------------------------------------*/

#include "varRhoIncompressibleTurbulenceModel.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(varRhoIncompressibleTurbulenceModel, 0);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::varRhoIncompressibleTurbulenceModel::varRhoIncompressibleTurbulenceModel
(
    const volScalarField& rho,
    const volVectorField& U,
    const surfaceScalarField& alphaRhoPhi,
    const surfaceScalarField& phi,
    const word& propertiesName
)
:
    turbulenceModel
    (
        U,
        alphaRhoPhi,
        phi,
        propertiesName
    ),
    rho_(rho)
{}

 Foam::tmp<Foam::surfaceScalarField>
 Foam::varRhoIncompressibleTurbulenceModel::phi() const
 {
     if (phi_.dimensions() == dimensionSet(0, 3, -1, 0, 0))
     {
         return phi_;
     }
     else
     {
         return phi_/fvc::interpolate(rho_);
     }
 }

Foam::tmp<Foam::volScalarField>
Foam::varRhoIncompressibleTurbulenceModel::mu() const
{
    return rho_*nu();
}


Foam::tmp<Foam::scalarField>
Foam::varRhoIncompressibleTurbulenceModel::mu(const label patchi) const
{
    return rho_*nu(patchi);
}


Foam::tmp<Foam::volScalarField>
Foam::varRhoIncompressibleTurbulenceModel::mut() const
{
    return rho_*nut();
}


Foam::tmp<Foam::scalarField>
Foam::varRhoIncompressibleTurbulenceModel::mut(const label patchi) const
{
    return rho_*nut(patchi);
}


Foam::tmp<Foam::volScalarField>
Foam::varRhoIncompressibleTurbulenceModel::muEff() const
{
    return rho_*nuEff();
}


Foam::tmp<Foam::scalarField>
Foam::varRhoIncompressibleTurbulenceModel::muEff(const label patchi) const
{
    return rho_*nuEff(patchi);
}


// ************************************************************************* //
