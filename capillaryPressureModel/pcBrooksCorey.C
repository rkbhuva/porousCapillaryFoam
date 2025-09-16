/*---------------------------------------------------------------------------*\
  =========                 |\
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |\
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |\
-------------------------------------------------------------------------------
    Copyright (C) 2024 Custom Implementation
\*---------------------------------------------------------------------------*/

#include "pcBrooksCorey.H"
#include "fvMesh.H"
#include "fvc.H"
#include "mathematicalConstants.H"
#include "addToRunTimeSelectionTable.H"
#include "IOstream.H"

namespace Foam
{
    // Define the run-time selection entry
    defineTypeNameAndDebug(pcBrooksCorey, 0);
    addToRunTimeSelectionTable(capillaryPressureModelBase, pcBrooksCorey, dictionary);
}

// capillary pressure calculation - Brooks Corey
void Foam::pcBrooksCorey::calculateCapillaryPressure()
{
        const scalar Srw = Srw_.value();
        const scalar lambda = max(lambda_.value(), SMALL);
        const scalar Pd = Pd_.value();
        const scalar maxPcLimit = maxPc_.value();

    forAll(pCapillary_, cellI)
    {
        scalar alphaW = (wettingPhase_ == 1) ? alpha1_[cellI] : alpha2_[cellI];
        scalar Se = (alphaW - Srw) / (1.0 - Srw);
        Se = max(1e-4, min(1.0 - 1e-4, Se));

        scalar pc = Pd * pow(Se, -1.0/lambda);
        pCapillary_[cellI] = max(min(pc, maxPcLimit), 0);

    }
    pCapillary_.correctBoundaryConditions();
}

Foam::tmp<Foam::volVectorField> Foam::pcBrooksCorey::capillaryPressureGradient() const
{
    tmp<volVectorField> tgradPCap
    (
        volVectorField::New
        (
            "gradPCap",
            mesh_,
            dimensionedVector("zero", dimPressure/dimLength, Zero)
        )
    );
    
    if (!enabled_) return tgradPCap;
    
    volVectorField& gradPCap = tgradPCap.ref();

    // Get the wetting phase field and its gradient
    const volScalarField& alphaW = (wettingPhase_ == 1) ? alpha1_ : alpha2_;
    tmp<volVectorField> gradAlpha = fvc::grad(alphaW);
    
    const scalar Srw = Srw_.value();
    const scalar lambda = max(lambda_.value(), SMALL);
    const scalar Pd = Pd_.value();
    const scalar maxPcLimit = maxPc_.value();
    
    forAll(gradPCap, cellI)
    {      
        scalar Se = (alphaW[cellI] - Srw) / (1.0 - Srw);
        Se = max(1e-4, min(1.0 - 1e-4, Se));

        scalar alphads = 1.0 / (1.0 - Srw); // (dSe/dα)
        vector galpha = gradAlpha()[cellI];

        if (pCapillary_[cellI] < maxPcLimit)
            gradPCap[cellI] = (-Pd * alphads / lambda * pow(Se, (-1.0/lambda) - 1.0)) * galpha; // ∇Pc = (dPc/dSe) × (dSe/dα) × ∇α
        else
            gradPCap[cellI] = vector::zero;
    }
    
    gradPCap.correctBoundaryConditions();
    return tgradPCap;
}

Foam::pcBrooksCorey::pcBrooksCorey
(
    const fvMesh& mesh,
    const volScalarField& alpha1,
    const volScalarField& alpha2,
    const dictionary& dict
)
:
    capillaryPressureModelBase(mesh, alpha1, alpha2, dict), // constructor
    Srw_("Srw", dimless, dict.getOrDefault<scalar>("Srw", 0.05)),
    lambda_("lambda", dimless, dict.getOrDefault<scalar>("lambda", 2.0)),
    Pd_("Pd", dimPressure, dict.getOrDefault<scalar>("Pd", 1.0))
{
    if (enabled_)
    {
        Info<< "Capillary pressure model enabled (Brooks-Corey)" << endl;
        Info<< "    Residual Saturation (Wetting): " << Srw_.value() << endl;
        Info<< "    Lambda: " << lambda_.value() << endl;
        Info<< "    Displacement Pressure: " << Pd_.value() << " Pa" << endl;
        Info<< "    Max Capillary Pressure: " << maxPc_.value() << " Pa" << endl;
    }

    calculateCapillaryPressure(); // Initial calculation
}

bool Foam::pcBrooksCorey::writeData(Ostream& os) const
{
    capillaryPressureModelBase::writeData(os);
    return os.good();
}

bool Foam::pcBrooksCorey::read(const dictionary& dict)
{
    capillaryPressureModelBase::read(dict);
    Srw_.value() = dict.getOrDefault<scalar>("Srw", Srw_.value());
    lambda_.value() = dict.getOrDefault<scalar>("lambda", lambda_.value());
    Pd_.value() = dict.getOrDefault<scalar>("Pd", Pd_.value());
    return true;
}

Foam::Ostream& operator<<(Foam::Ostream& os, const Foam::pcBrooksCorey& model)
{
    model.writeData(os);
    return os;
}