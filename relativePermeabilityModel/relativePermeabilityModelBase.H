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

#ifndef relativePermeabilityModelBase_H
#define relativePermeabilityModelBase_H

#include "fvMesh.H"
#include "volFields.H"
#include "dictionary.H"
#include "dimensionedScalar.H"
#include "Switch.H"
#include "IOobject.H"
#include "autoPtr.H"
#include "className.H"
#include "runTimeSelectionTables.H"

namespace Foam
{

class relativePermeabilityModelBase
{
protected:
    const fvMesh& mesh_;
    const volScalarField& alpha1_;
    const volScalarField& alpha2_;
    dictionary coeffDict_;

    volScalarField Kr1_;
    volScalarField Kr2_;
    Switch enabled_;
    label wettingPhase_;

    dimensionedScalar D_;

private:
    relativePermeabilityModelBase(const relativePermeabilityModelBase&) = delete;
    void operator=(const relativePermeabilityModelBase&) = delete;

public:

    TypeName("relativePermeabilityModelBase");

    declareRunTimeSelectionTable
    (
        autoPtr,
        relativePermeabilityModelBase,
        dictionary,
        (
            const fvMesh& mesh,
            const volScalarField& alpha1,
            const volScalarField& alpha2,
            const dictionary& dict
        ),
        (mesh, alpha1, alpha2, dict)
    );

    relativePermeabilityModelBase
    (
        const fvMesh& mesh,
        const volScalarField& alpha1,
        const volScalarField& alpha2,
        const dictionary& dict
    );

    virtual ~relativePermeabilityModelBase() = default;

    static autoPtr<relativePermeabilityModelBase> New
    (
        const fvMesh& mesh,
        const volScalarField& alpha1,
        const volScalarField& alpha2,
        const dictionary& dict
    );

    const volScalarField& Kr1() const { return Kr1_; }
    const volScalarField& Kr2() const { return Kr2_; }
    bool enabled() const { return enabled_; }
    label wettingPhase() const { return wettingPhase_; }
    const dimensionedScalar& D() const { return D_; }

    virtual void calculateRelativePermeability() = 0;
    virtual void update();

    virtual bool writeData(Ostream& os) const;
    virtual bool read(const dictionary& dict);
};

Ostream& operator<<(Ostream&, const relativePermeabilityModelBase&);

}

#endif
