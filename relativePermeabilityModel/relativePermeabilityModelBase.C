/*---------------------------------------------------------------------------*\
  =========                 |\
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |\
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |\
-------------------------------------------------------------------------------
    Copyright (C) 2024 Custom Implementation
\*---------------------------------------------------------------------------*/

#include "relativePermeabilityModelBase.H"
#include "addToRunTimeSelectionTable.H"
#include "fvc.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "IOstreams.H"

namespace Foam
{
    defineTypeNameAndDebug(relativePermeabilityModelBase, 0);
    defineRunTimeSelectionTable(relativePermeabilityModelBase, dictionary);
}

Foam::relativePermeabilityModelBase::relativePermeabilityModelBase
(
    const fvMesh& mesh,
    const volScalarField& alpha1,
    const volScalarField& alpha2,
    const dictionary& dict
)
:
    mesh_(mesh),
    alpha1_(alpha1),
    alpha2_(alpha2),
    coeffDict_(dict),

    Kr1_
    (
        IOobject
        (
            "Kr1",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("zero", dimless, 0.0)
    ),

    Kr2_
    (
        IOobject
        (
            "Kr2",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("zero", dimless, 0.0)
    ),

    enabled_(dict.getOrDefault<Switch>("enabled", true)),
    wettingPhase_(dict.getOrDefault<label>("wettingPhase", 1)),

    D_("D", dimensionSet(0, -2, 0, 0, 0, 0, 0), dict.getOrDefault<scalar>("D", 1e12))
{}


Foam::autoPtr<Foam::relativePermeabilityModelBase>
Foam::relativePermeabilityModelBase::New
(
    const fvMesh& mesh,
    const volScalarField& alpha1,
    const volScalarField& alpha2,
    const dictionary& globalDict
)
{
    const dictionary& dict(globalDict.subDict("relativePermeabilityModel"));

    word modelType(dict.lookup("type"));

    Info<< "Selecting relative permeability model " << modelType << endl;

    auto* ctorPtr = dictionaryConstructorTable(modelType);

    if (!ctorPtr)
    {
        FatalErrorInFunction
            << "Unknown relativePermeabilityModel type " << modelType << nl
            << "Valid relativePermeabilityModel types are:" << nl
            << dictionaryConstructorTablePtr_->toc()
            << exit(FatalError);
    }

    return autoPtr<relativePermeabilityModelBase>(ctorPtr(mesh, alpha1, alpha2, dict));
}


void Foam::relativePermeabilityModelBase::update()
{
    if (enabled_)
    {
        calculateRelativePermeability();
    }
}

bool Foam::relativePermeabilityModelBase::writeData(Ostream& os) const
{
    Kr1_.write();
    Kr2_.write();
    return os.good();
}

bool Foam::relativePermeabilityModelBase::read(const dictionary& dict)
{
    enabled_ = dict.getOrDefault<Switch>("enabled", enabled_);
    wettingPhase_ = dict.getOrDefault<label>("wettingPhase", wettingPhase_);
    D_.value() = dict.getOrDefault<scalar>("D", D_.value());
    return true;
}

Foam::Ostream& Foam::operator<<(Foam::Ostream& os, const Foam::relativePermeabilityModelBase& model)
{
    model.writeData(os);
    return os;
}
