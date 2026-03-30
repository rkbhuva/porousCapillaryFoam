/*---------------------------------------------------------------------------*\
  =========                 |\
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |\
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |\
-------------------------------------------------------------------------------
    Copyright (C) 2024 Custom Implementation
\*---------------------------------------------------------------------------*/

#include "krCorey.H"
#include "addToRunTimeSelectionTable.H"
#include "fvc.H"
#include "IOstream.H"

namespace Foam
{
    defineTypeNameAndDebug(krCorey, 0);
    addToRunTimeSelectionTable(relativePermeabilityModelBase, krCorey, dictionary);
}

void Foam::krCorey::calculateRelativePermeability()
{
    forAll(Kr1_, cellI)
    {
        scalar alphaW = (wettingPhase_ == 1) ? alpha1_[cellI] : alpha2_[cellI];

        scalar Sr = Sr_.value();
        scalar n = n_.value();

        scalar Se = (alphaW - Sr) / (1.0 - Sr);
        Se = max(min(Se, 1.0 - 1e-6), 1e-6);

        scalar Krw_val, Krnw_val;

        Krw_val = pow(Se, n);
        Krnw_val = pow((1 - Se), n);
  
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

Foam::krCorey::krCorey
(
    const fvMesh& mesh,
    const volScalarField& alpha1,
    const volScalarField& alpha2,
    const dictionary& dict
)
:
    relativePermeabilityModelBase(mesh, alpha1, alpha2, dict),
    Sr_("Sr", dimless, dict.getOrDefault<scalar>("Sr", 0.05)),
    n_("n", dimless, dict.getOrDefault<scalar>("n", 3.0))
{
    if (enabled_)
    {
        Info<< "Relative permeability model enabled (Corey)" << endl;
        Info<< "    Residual Saturation (Wetting): " << Sr_.value() << endl;
        Info<< "    n: " << n_.value() << endl;
        Info<< "    Darcy coefficient: " << D_.value() << " 1/m^2" << endl;
    }
    calculateRelativePermeability();
}

bool Foam::krCorey::writeData(Ostream& os) const
{
    relativePermeabilityModelBase::writeData(os);
    return os.good();
}

bool Foam::krCorey::read(const dictionary& dict)
{
    relativePermeabilityModelBase::read(dict);
    Sr_.value() = dict.getOrDefault("Sr", Sr_.value());
    n_.value() = dict.getOrDefault("n", n_.value());
    return true;
}

Foam::Ostream& Foam::operator<<(Foam::Ostream& os, const krCorey& model)
{
    model.writeData(os);
    return os;
}
