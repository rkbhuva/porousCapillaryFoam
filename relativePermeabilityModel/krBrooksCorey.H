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

#ifndef krbrooksCorey_H
#define krbrooksCorey_H

#include "relativePermeabilityModelBase.H"

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class krBrooksCorey Declaration
\*---------------------------------------------------------------------------*/

class krBrooksCorey
:
    public relativePermeabilityModelBase
{
    dimensionedScalar Srw_;
    dimensionedScalar lambda_;

    void calculateRelativePermeability() override;

    krBrooksCorey(const krBrooksCorey&) = delete;
    void operator=(const krBrooksCorey&) = delete;

public:

    TypeName("krBrooksCorey");

    krBrooksCorey
    (
        const fvMesh& mesh,
        const volScalarField& alpha1,
        const volScalarField& alpha2,
        const dictionary& dict
    );

    virtual ~krBrooksCorey() = default;

    bool writeData(Ostream& os) const override;
    bool read(const dictionary& dict) override;
};

Foam::Ostream& operator<<(Foam::Ostream&, const krBrooksCorey&);

}

#endif
