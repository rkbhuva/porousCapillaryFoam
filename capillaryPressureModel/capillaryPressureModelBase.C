/*---------------------------------------------------------------------------*\
  =========                 |\
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |\
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |\
-------------------------------------------------------------------------------
    Copyright (C) 2024 Custom Implementation
\*---------------------------------------------------------------------------*/

#include "capillaryPressureModelBase.H"
#include "fvMesh.H"
#include "fvc.H"
#include "addToRunTimeSelectionTable.H"
#include "IOstream.H"

namespace Foam
{
    defineTypeNameAndDebug(capillaryPressureModelBase, 0);
    defineRunTimeSelectionTable(capillaryPressureModelBase, dictionary);
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::capillaryPressureModelBase::capillaryPressureModelBase
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
    pCapillary_
    (
        IOobject
        (
            "pCapillary",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar("pCapillary", dimPressure, 0.0)
    ),
    enabled_(dict.getOrDefault<Switch>("enabled", true)),
    wettingPhase_(dict.getOrDefault<label>("wettingPhase", 1)),
    maxPc_
    (
        "maxPc", dimPressure,
        dict.getOrDefault<scalar>("maxPc", 1e5)
    )
    
{}


// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::capillaryPressureModelBase> Foam::capillaryPressureModelBase::New
(
    const fvMesh& mesh,
    const volScalarField& alpha1,
    const volScalarField& alpha2,
    const dictionary& globalDict
)
{
    const dictionary& dict(globalDict.subDict("capillaryPressureModel"));

    word modelType(dict.lookup("type"));

    Info<< "Selecting capillary pressure model " << modelType << endl;

    auto* ctorPtr = dictionaryConstructorTable(modelType);

    if (!ctorPtr)
    {
        FatalErrorInFunction
            << "Unknown capillaryPressureModel type " << modelType << nl
            << "Valid capillaryPressureModel types are:" << nl
            << dictionaryConstructorTablePtr_->toc()
            << exit(FatalError);
    }

    return autoPtr<capillaryPressureModelBase>(ctorPtr(mesh, alpha1, alpha2, dict));
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //
void Foam::capillaryPressureModelBase::update()
{
    if (!enabled_) return;

    calculateCapillaryPressure();
}

Foam::tmp<Foam::volVectorField>
Foam::capillaryPressureModelBase::capillaryPressureMomentum() const
{
    return capillaryPressureGradient();
}

bool Foam::capillaryPressureModelBase::writeData(Ostream& os) const
{
    pCapillary_.write();
    return os.good();
}

bool Foam::capillaryPressureModelBase::read(const dictionary& dict)
{
    enabled_ = dict.getOrDefault<Switch>("enabled", enabled_);
    wettingPhase_ = dict.getOrDefault<label>("wettingPhase", wettingPhase_);
    maxPc_.value() = dict.getOrDefault<scalar>("maxPc", maxPc_.value());
    return true;
}

Foam::Ostream& operator<<(Foam::Ostream& os, const Foam::capillaryPressureModelBase& model)
{
    model.writeData(os);
    return os;
}