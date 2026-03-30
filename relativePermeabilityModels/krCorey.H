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

#ifndef krCorey_H
#define krCorey_H

#include "relativePermeabilityModelBase.H"

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class krCorey Declaration
\*---------------------------------------------------------------------------*/

class krCorey
:
    public relativePermeabilityModelBase
{
    dimensionedScalar Sr_;
    dimensionedScalar n_;

    void calculateRelativePermeability() override;

    krCorey(const krCorey&) = delete;
    void operator=(const krCorey&) = delete;

public:

    TypeName("krCorey");

    krCorey
    (
        const fvMesh& mesh,
        const volScalarField& alpha1,
        const volScalarField& alpha2,
        const dictionary& dict
    );

    virtual ~krCorey() = default;

    bool writeData(Ostream& os) const override;
    bool read(const dictionary& dict) override;
};

Foam::Ostream& operator<<(Foam::Ostream&, const krCorey&);

}

#endif
