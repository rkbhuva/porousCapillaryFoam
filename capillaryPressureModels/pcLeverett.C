/*---------------------------------------------------------------------------*\
  =========                 |\
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |\
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |\
-------------------------------------------------------------------------------
    Copyright (C) 2024 Custom Implementation
\*---------------------------------------------------------------------------*/

#include "pcLeverett.H"
#include "fvMesh.H"
#include "fvc.H"
#include "mathematicalConstants.H"
#include "addToRunTimeSelectionTable.H"
#include "IOstream.H"

namespace Foam
{
    // Define the run-time selection entry
    defineTypeNameAndDebug(pcLeverett, 0);
    addToRunTimeSelectionTable(capillaryPressureModelBase, pcLeverett, dictionary);
}

Foam::tmp<Foam::volVectorField> Foam::pcLeverett::capillaryPressureGradient() const
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
    
    const scalar sigma = sigma_.value();
    const scalar thetaDeg = theta_.value();
    const scalar thetaRad = thetaDeg * Foam::constant::mathematical::pi / 180.0;
    const scalar D = D_.value();
    const scalar phi = phi_.value();
    const scalar maxPcLimit = maxPc_.value();
    
    forAll(gradPCap, cellI)
    {  
        scalar alphaSe = alphaW[cellI];
        
        // 1. Calculate Pc
        scalar pcVal = (sigma * cos(thetaRad)) * sqrt(D * phi) 
            * (thetaDeg < 90.0 ? (1.417 * (1.0 - alphaSe)) - (2.120 * pow((1.0 - alphaSe), 2)) + (1.263 * pow((1.0 - alphaSe), 3))
            : (1.417 * alphaSe) - (2.120 * pow(alphaSe, 2)) + (1.263 * pow(alphaSe, 3)) );
            
        pc[cellI] = max(min(pcVal, maxPcLimit), 0.0);

        // 2. Calculate Gradient
        if (pc[cellI] < maxPcLimit)
        {
            vector galpha = gradAlpha()[cellI];
            scalar dpcdalpha = (sigma * cos(thetaRad)) * sqrt(D * phi) 
                * (thetaDeg < 90.0 ? (-1.417) + (2.0 * 2.120 * (1.0 - alphaSe)) - (3.0 * 1.263 * pow((1.0 - alphaSe), 2))
                : 1.417 - (2.0 * 2.120 * alphaSe) + (3.0 * 1.263 * pow(alphaSe, 2)) );
                
            gradPCap[cellI] = dpcdalpha * galpha; 
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

Foam::pcLeverett::pcLeverett
(
    const fvMesh& mesh,
    const volScalarField& alpha1,
    const volScalarField& alpha2,
    const dictionary& dict
)
:
    capillaryPressureModelBase(mesh, alpha1, alpha2, dict), // constructor
    theta_("theta",dimless,dict.getOrDefault<scalar>("theta", 60.0)),
    sigma_("sigma",dimForce/dimLength,dict.getOrDefault<scalar>("sigma", 0.072)),
    phi_("phi",dimless,dict.getOrDefault<scalar>("phi", 0.3)),
    D_("D", dimensionSet(0, -2, 0, 0, 0, 0, 0), dict.getOrDefault<scalar>("D", 1e12))
{
    if (enabled_)
    {
        Info<< "Capillary pressure model enabled (Leverett):" << endl;
        Info<< "    Surface tension: " << sigma_.value() << " N/m" << endl;
        Info<< "    Contact Angle: " << theta_.value() << " deg" << endl;
        Info<< "    Porosity: " << phi_.value() << endl;
        Info<< "    Darcy coefficient: " << D_.value() << " 1/m^2" << endl;
        Info<< "    Max Capillary Pressure: " << maxPc_.value() << " Pa" << endl; 
    }

    (void)capillaryPressureGradient(); // Initial calculation
}

bool Foam::pcLeverett::writeData(Ostream& os) const
{
    capillaryPressureModelBase::writeData(os);
    return os.good();
}

bool Foam::pcLeverett::read(const dictionary& dict)
{
    capillaryPressureModelBase::read(dict);
    sigma_.value() = dict.getOrDefault<scalar>("sigma", sigma_.value());
    phi_.value() = dict.getOrDefault<scalar>("phi", phi_.value());
    theta_.value() = dict.getOrDefault<scalar>("theta", theta_.value());
    D_.value() = dict.getOrDefault<scalar>("D", D_.value());
    return true;
}

Foam::Ostream& operator<<(Foam::Ostream& os, const Foam::pcLeverett& model)
{
    model.writeData(os);
    return os;
}
