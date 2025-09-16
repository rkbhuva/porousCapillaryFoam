/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2024 Custom Implementation
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

\*---------------------------------------------------------------------------*/

#ifndef pcVanGenuchten_H
#define pcVanGenuchten_H

#include "capillaryPressureModelBase.H"

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class pcVanGenuchten Declaration
\*---------------------------------------------------------------------------*/

class pcVanGenuchten
:
    public capillaryPressureModelBase
{
    // Private Data (specific to pcVanGenuchten)
    dimensionedScalar Srw_;
    dimensionedScalar m_;
    dimensionedScalar pc0_;

    // Private Member Functions
    void calculateCapillaryPressure() override;

    //- No copy construct
    pcVanGenuchten(const pcVanGenuchten&) = delete;

    //- No copy assignment
    void operator=(const pcVanGenuchten&) = delete;

public:

    TypeName("pcVanGenuchten");

    // Constructor
    pcVanGenuchten
    (
        const fvMesh& mesh,
        const volScalarField& alpha1,
        const volScalarField& alpha2,
        const dictionary& dict
    );

    virtual ~pcVanGenuchten() = default;
    virtual tmp<volVectorField> capillaryPressureGradient() const override;

    // I-O
    bool writeData(Ostream& os) const override;
    bool read(const dictionary& dict) override;
};

// Ostream operator (for Info<< model)
Foam::Ostream& operator<<(Foam::Ostream&, const pcVanGenuchten&);
}

#endif