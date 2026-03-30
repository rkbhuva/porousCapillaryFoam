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
    volScalarField& pc = const_cast<volScalarField&>(pc_);

    // Get the wetting phase field and its gradient
    const volScalarField& alphaW = (wettingPhase_ == 1) ? alpha1_ : alpha2_;
    tmp<volVectorField> gradAlpha = fvc::grad(alphaW);
    
    const scalar Sr = Sr_.value();
    const scalar lambda = max(lambda_.value(), SMALL);
    const scalar pe = pe_.value();
    const scalar maxPcLimit = maxPc_.value();
    
    forAll(gradPCap, cellI)
    {
        scalar alphaCell = alphaW[cellI];      
        scalar Se = (alphaCell - Sr) / (1.0 - Sr);
        Se = max(min(Se, 1.0 - 1e-6), 1e-6);

        scalar pcVal = pe * pow(Se, -1.0/lambda);
        pc[cellI] = min(pcVal, maxPcLimit);

        if (pc[cellI] < maxPcLimit)
        {
            scalar alphads = 1.0 / (1.0 - Sr); // (dSe/dα)
            vector galpha = gradAlpha()[cellI];
            gradPCap[cellI] = (-pe * alphads / lambda * pow(Se, (-1.0/lambda) - 1.0)) * galpha; // ∇Pc = (dPc/dSe) × (dSe/dα) × ∇α
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

Foam::pcBrooksCorey::pcBrooksCorey
(
    const fvMesh& mesh,
    const volScalarField& alpha1,
    const volScalarField& alpha2,
    const dictionary& dict
)
:
    capillaryPressureModelBase(mesh, alpha1, alpha2, dict), // constructor
    Sr_("Sr", dimless, dict.getOrDefault<scalar>("Sr", 0.05)),
    lambda_("lambda", dimless, dict.getOrDefault<scalar>("lambda", 2.0)),
    pe_("pe", dimPressure, dict.getOrDefault<scalar>("pe", 1.0))
{
    if (enabled_)
    {
        Info<< "Capillary pressure model enabled (Brooks-Corey)" << endl;
        Info<< "    Residual Saturation (Wetting): " << Sr_.value() << endl;
        Info<< "    Lambda: " << lambda_.value() << endl;
        Info<< "    Entry Pressure: " << pe_.value() << " Pa" << endl;
        Info<< "    Max Capillary Pressure: " << maxPc_.value() << " Pa" << endl;
    }

    (void)capillaryPressureGradient(); // Initial calculation
}

bool Foam::pcBrooksCorey::writeData(Ostream& os) const
{
    capillaryPressureModelBase::writeData(os);
    return os.good();
}

bool Foam::pcBrooksCorey::read(const dictionary& dict)
{
    capillaryPressureModelBase::read(dict);
    Sr_.value() = dict.getOrDefault<scalar>("Sr", Sr_.value());
    lambda_.value() = dict.getOrDefault<scalar>("lambda", lambda_.value());
    pe_.value() = dict.getOrDefault<scalar>("pe", pe_.value());
    return true;
}

Foam::Ostream& operator<<(Foam::Ostream& os, const Foam::pcBrooksCorey& model)
{
    model.writeData(os);
    return os;
}
