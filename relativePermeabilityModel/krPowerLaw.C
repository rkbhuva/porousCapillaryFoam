/*---------------------------------------------------------------------------*\
  =========                 |\
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |\
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |\
-------------------------------------------------------------------------------
    Copyright (C) 2024 Custom Implementation
\*---------------------------------------------------------------------------*/

#include "krPowerLaw.H"
#include "addToRunTimeSelectionTable.H"
#include "fvc.H"
#include "IOstream.H"

namespace Foam
{
    defineTypeNameAndDebug(krPowerLaw, 0);
    addToRunTimeSelectionTable(relativePermeabilityModelBase, krPowerLaw, dictionary);
}

void Foam::krPowerLaw::calculateRelativePermeability()
{
    forAll(Kr1_, cellI)
    {
        scalar alphaW = (wettingPhase_ == 1) ? alpha1_[cellI] : alpha2_[cellI];

        scalar pwr = pwr_.value();
        scalar Krw_val, Krnw_val;

        Krw_val = pow(alphaW, pwr);
        Krnw_val = pow((1 - alphaW), pwr);
  
        Krw_val = max(Krw_val, 1e-6);
        Krnw_val = max(Krnw_val, 1e-6);

        if (wettingPhase_ == 1)
        {
            Kr1_[cellI] = Krw_val;
            Kr2_[cellI] = Krnw_val;
        }
        else
        {
            Kr1_[cellI] = Krnw_val;
            Kr2_[cellI] = Krw_val;
        }
    }
    Kr1_.correctBoundaryConditions();
    Kr2_.correctBoundaryConditions();
}

Foam::krPowerLaw::krPowerLaw
(
    const fvMesh& mesh,
    const volScalarField& alpha1,
    const volScalarField& alpha2,
    const dictionary& dict
)
:
    relativePermeabilityModelBase(mesh, alpha1, alpha2, dict),
    pwr_("pwr", dimless, dict.getOrDefault<scalar>("pwr", 3.0))
{
    if (enabled_)
    {
        Info<< "Relative permeability model enabled (Power Law)" << endl;
        Info<< "    Power: " << pwr_.value() << endl;
        Info<< "    Darcy coefficient: " << D_.value() << " 1/m^2" << endl;
    }
    calculateRelativePermeability();
}

bool Foam::krPowerLaw::writeData(Ostream& os) const
{
    relativePermeabilityModelBase::writeData(os);
    return os.good();
}

bool Foam::krPowerLaw::read(const dictionary& dict)
{
    relativePermeabilityModelBase::read(dict);
    pwr_.value() = dict.getOrDefault("pwr", pwr_.value());
    return true;
}

Foam::Ostream& Foam::operator<<(Foam::Ostream& os, const krPowerLaw& model)
{
    model.writeData(os);
    return os;
}
