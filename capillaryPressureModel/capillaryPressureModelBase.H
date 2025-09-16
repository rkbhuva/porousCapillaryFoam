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

#ifndef capillaryPressureModelBase_H
#define capillaryPressureModelBase_H

#include "fvMesh.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "dictionary.H"
#include "dimensionedScalar.H"
#include "dimensionedVector.H"
#include "Switch.H"
#include "IOobject.H"
#include "autoPtr.H"
#include "className.H"
#include "runTimeSelectionTables.H"

namespace Foam
{

class capillaryPressureModelBase
{
protected:
    // Protected Data (common to all models)
    const fvMesh& mesh_;
    const volScalarField& alpha1_;
    const volScalarField& alpha2_;
    dictionary coeffDict_;

    volScalarField pCapillary_;
    Switch enabled_;
    label wettingPhase_;

    dimensionedScalar maxPc_;

private:
    //- No copy construct
    capillaryPressureModelBase(const capillaryPressureModelBase&) = delete;

    //- No copy assignment
    void operator=(const capillaryPressureModelBase&) = delete;

public:

    TypeName("capillaryPressureModelBase");

    declareRunTimeSelectionTable
    (
        autoPtr,
        capillaryPressureModelBase,
        dictionary,
        (
            const fvMesh& mesh,
            const volScalarField& alpha1,
            const volScalarField& alpha2,
            const dictionary& dict
        ),
        (mesh, alpha1, alpha2, dict)
    );

    // Constructor
    capillaryPressureModelBase
    (
        const fvMesh& mesh,
        const volScalarField& alpha1,
        const volScalarField& alpha2,
        const dictionary& dict
    );

    virtual ~capillaryPressureModelBase() = default;

    // Selector
    static autoPtr<capillaryPressureModelBase> New
    (
        const fvMesh& mesh,
        const volScalarField& alpha1,
        const volScalarField& alpha2,
        const dictionary& dict
    );

    // Accessors (common to all models)
    const volScalarField& pCapillary() const { return pCapillary_; }
    bool enabled() const { return enabled_; }
    label wettingPhase() const { return wettingPhase_; }
    const dimensionedScalar& maxPc() const { return maxPc_; }


    virtual void calculateCapillaryPressure() = 0;
    virtual tmp<volVectorField> capillaryPressureGradient() const = 0;
    virtual void update();

    virtual tmp<volVectorField> capillaryPressureMomentum() const;
    virtual bool writeData(Ostream& os) const;
    virtual bool read(const dictionary& dict);
};

Foam::Ostream& operator<<(Foam::Ostream&, const capillaryPressureModelBase&);

}

#endif