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
    volScalarField& pc = const_cast<volScalarField&>(pc_);

    // Get the wetting phase field and its gradient
    const volScalarField& alphaW = (wettingPhase_ == 1) ? alpha1_ : alpha2_;
    tmp<volVectorField> gradAlpha = fvc::grad(alphaW);
    
    const scalar Sr = Sr_.value();
    const scalar m = max(m_.value(), SMALL);
    const scalar p0 = p0_.value();
    const scalar maxPcLimit = maxPc_.value();
    
    forAll(gradPCap, cellI)
    {
        scalar alphaCell = alphaW[cellI];      
        scalar Se = (alphaCell - Sr) / (1.0 - Sr);
        Se = max(min(Se, 1.0 - 1e-6), 1e-6);

        scalar pcVal = p0 * pow(((pow(Se, -1.0/m)) - 1.0), 1-m);
        pc[cellI] = min(pcVal, maxPcLimit);

        if (pc[cellI] < maxPcLimit)
        {
            scalar alphads = 1.0 / (1.0 - Sr); // (dSe/dα)
            vector galpha = gradAlpha()[cellI];
            gradPCap[cellI] = (p0 * alphads * (m - 1.0)/m * pow(Se, -(1.0 + m)/m) * pow((pow(Se, -1.0/m) - 1.0), -m)) * galpha; 
        }
        else
        {
            gradPCap[cellI] = vector::zero;
        }
    }
    
    pc.correctBoundaryConditions();
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
    Sr_("Sr", dimless, dict.getOrDefault<scalar>("Sr", 0.05)),
    m_("m", dimless, dict.getOrDefault<scalar>("m", 0.5)),
    p0_("p0", dimPressure, dict.getOrDefault<scalar>("p0", 1.0))
{
    if (enabled_)
    {
        Info<< "Capillary pressure model enabled (Van Genuchten)" << endl;
        Info<< "    Residual Saturation (Wetting): " << Sr_.value() << endl;
        Info<< "    m: " << m_.value() << endl;
        Info<< "    Entry Pressure: " << p0_.value() << " Pa" << endl;
        Info<< "    Max Capillary Pressure: " << maxPc_.value() << " Pa" << endl;
    }

    (void)capillaryPressureGradient(); // Initial calculation
}

bool Foam::pcVanGenuchten::writeData(Ostream& os) const
{
    capillaryPressureModelBase::writeData(os);
    return os.good();
}

bool Foam::pcVanGenuchten::read(const dictionary& dict)
{
    capillaryPressureModelBase::read(dict);
    Sr_.value() = dict.getOrDefault<scalar>("Sr", Sr_.value());
    m_.value() = dict.getOrDefault<scalar>("m", m_.value());
    p0_.value() = dict.getOrDefault<scalar>("p0", p0_.value());
    return true;
}

Foam::Ostream& operator<<(Foam::Ostream& os, const Foam::pcVanGenuchten& model)
{
    model.writeData(os);
    return os;
}
