/*---------------------------------------------------------------------------*\
  =========                 |\
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |\
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |\
-------------------------------------------------------------------------------
    Copyright (C) 2024 Custom Implementation
\*---------------------------------------------------------------------------*/

#include "pcVanGenuchten.H"
#include "fvMesh.H"
#include "fvc.H"
#include "mathematicalConstants.H"
#include "addToRunTimeSelectionTable.H"
#include "IOstream.H"

namespace Foam
{
    // Define the run-time selection entry
    defineTypeNameAndDebug(pcVanGenuchten, 0);
    addToRunTimeSelectionTable(capillaryPressureModelBase, pcVanGenuchten, dictionary);
}

// capillary pressure calculation - pcVanGenuchten
void Foam::pcVanGenuchten::calculateCapillaryPressure()
{
    const scalar Srw = Srw_.value();
    const scalar m = max(m_.value(), SMALL);
    const scalar p0 = pc0_.value();
    const scalar maxPcLimit = maxPc_.value();

    forAll(pCapillary_, cellI)
    {
        scalar alphaW = (wettingPhase_ == 1) ? alpha1_[cellI] : alpha2_[cellI];
        scalar Se = (alphaW - Srw) / (1.0 - Srw);
        Se = max(1e-4, min(1.0 - 1e-4, Se));
        
        scalar pc = p0 * pow(((pow(Se, -1.0/m)) - 1.0), 1-m);
        pCapillary_[cellI] = max(min(pc, maxPcLimit), 0.0);
    }
    pCapillary_.correctBoundaryConditions();
}

Foam::tmp<Foam::volVectorField> Foam::pcVanGenuchten::capillaryPressureGradient() const
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
    const scalar m = max(m_.value(), SMALL);
    const scalar p0 = pc0_.value();
    const scalar maxPcLimit = maxPc_.value();
    
    forAll(gradPCap, cellI)
    {      
        scalar Se = (alphaW[cellI] - Srw) / (1.0 - Srw);
        Se = max(1e-4, min(1.0 - 1e-4, Se));

        scalar alphads = 1.0 / (1.0 - Srw); // (dSe/dα)
        vector galpha = gradAlpha()[cellI];

        if (pCapillary_[cellI] < maxPcLimit)
            gradPCap[cellI] = (p0 * alphads * (m - 1.0)/m * pow(Se, -(1.0 + m)/m) * pow((pow(Se, -1.0/m) - 1.0), -m)) * galpha; // ∇Pc = (dPc/dSe) × (dSe/dα) × ∇α
        else
            gradPCap[cellI] = vector::zero;
    }
    
    gradPCap.correctBoundaryConditions();
    return tgradPCap;
}

Foam::pcVanGenuchten::pcVanGenuchten
(
    const fvMesh& mesh,
    const volScalarField& alpha1,
    const volScalarField& alpha2,
    const dictionary& dict
)
:
    capillaryPressureModelBase(mesh, alpha1, alpha2, dict), // constructor
    Srw_("Srw", dimless, dict.getOrDefault<scalar>("Srw", 0.05)),
    m_("m", dimless, dict.getOrDefault<scalar>("m", 0.5)),
    pc0_("pc0", dimPressure, dict.getOrDefault<scalar>("pc0", 1.0))
{
    if (enabled_)
    {
        Info<< "Capillary pressure model enabled (Van Genuchten)" << endl;
        Info<< "    Residual Saturation (Wetting): " << Srw_.value() << endl;
        Info<< "    m: " << m_.value() << endl;
        Info<< "    Pressure Pc0: " << pc0_.value() << " Pa" << endl;
        Info<< "    Max Capillary Pressure: " << maxPc_.value() << " Pa" << endl;
    }

    calculateCapillaryPressure(); // Initial calculation
}

bool Foam::pcVanGenuchten::writeData(Ostream& os) const
{
    capillaryPressureModelBase::writeData(os);
    return os.good();
}

bool Foam::pcVanGenuchten::read(const dictionary& dict)
{
    capillaryPressureModelBase::read(dict);
    Srw_.value() = dict.getOrDefault<scalar>("Srw", Srw_.value());
    m_.value() = dict.getOrDefault<scalar>("m", m_.value());
    pc0_.value() = dict.getOrDefault<scalar>("pc0", pc0_.value());
    return true;
}

Foam::Ostream& operator<<(Foam::Ostream& os, const Foam::pcVanGenuchten& model)
{
    model.writeData(os);
    return os;
}