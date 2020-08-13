/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2013-2016 OpenFOAM Foundation
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

#include "VarRhoIncompressibleTurbulenceModel.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class TransportModel>
Foam::VarRhoIncompressibleTurbulenceModel<TransportModel>::
VarRhoIncompressibleTurbulenceModel
(
    const word& type,
    const geometricOneField& alpha,
    const volScalarField& rho,
    const volVectorField& U,
    const surfaceScalarField& alphaRhoPhi,
    const surfaceScalarField& phi,
    const TransportModel& transport,
    const word& propertiesName
)
:
    TurbulenceModel
    <
        geometricOneField,
        volScalarField,
        varRhoIncompressibleTurbulenceModel,
        TransportModel
    >
    (
        alpha,
        rho,
        U,
        alphaRhoPhi,
        phi,
        transport,
        propertiesName
    )
{}


// * * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * //

template<class TransportModel>
Foam::autoPtr<Foam::VarRhoIncompressibleTurbulenceModel<TransportModel>>
Foam::VarRhoIncompressibleTurbulenceModel<TransportModel>::New
(
    const volScalarField& rho,
    const volVectorField& U,
    const surfaceScalarField& alphaRhoPhi,
    const surfaceScalarField& phi,
    const TransportModel& transport,
    const word& propertiesName
)
{
    return autoPtr<VarRhoIncompressibleTurbulenceModel>
    (
        static_cast<VarRhoIncompressibleTurbulenceModel*>(
        TurbulenceModel
        <
            geometricOneField,
            volScalarField,
            varRhoIncompressibleTurbulenceModel,
            TransportModel
        >::New
        (
            geometricOneField(),
            rho,
            U,
            alphaRhoPhi,
            phi,
            transport,
            propertiesName
        ).ptr())
    );
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class TransportModel>
Foam::tmp<Foam::volSymmTensorField>
Foam::VarRhoIncompressibleTurbulenceModel<TransportModel>::devReff() const
{
    return devRhoReff();
}


template<class TransportModel>
Foam::tmp<Foam::fvVectorMatrix>
Foam::VarRhoIncompressibleTurbulenceModel<TransportModel>::divDevReff
(
    volVectorField& U
) const
{
    return divDevRhoReff(U);
}


template<class TransportModel>
Foam::tmp<Foam::volSymmTensorField>
Foam::VarRhoIncompressibleTurbulenceModel<TransportModel>::
devRhoReff() const
{
    NotImplemented;

    return devReff();
}


template<class TransportModel>
Foam::tmp<Foam::fvVectorMatrix>
Foam::VarRhoIncompressibleTurbulenceModel<TransportModel>::
divDevRhoReff
(
    volVectorField& U
) const
{
    NotImplemented;

    return divDevReff(U);
}


template<class TransportModel>
Foam::tmp<Foam::fvVectorMatrix>
Foam::VarRhoIncompressibleTurbulenceModel<TransportModel>::
divDevRhoReff
(
    const volScalarField& rho,
    volVectorField& U
) const
{
    NotImplemented;

    return divDevReff(U);
}


// ************************************************************************* //
