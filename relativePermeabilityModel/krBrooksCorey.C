/*---------------------------------------------------------------------------*\
  =========                 |\
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |\
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |\
-------------------------------------------------------------------------------
    Copyright (C) 2024 Custom Implementation
\*---------------------------------------------------------------------------*/

#include "krBrooksCorey.H"
#include "addToRunTimeSelectionTable.H"
#include "fvc.H"
#include "IOstream.H"

namespace Foam
{
    defineTypeNameAndDebug(krBrooksCorey, 0);
    addToRunTimeSelectionTable(relativePermeabilityModelBase, krBrooksCorey, dictionary);
}

void Foam::krBrooksCorey::calculateRelativePermeability()
{
    forAll(Kr1_, cellI)
    {
        scalar alphaW = (wettingPhase_ == 1) ? alpha1_[cellI] : alpha2_[cellI];

        scalar Srw = Srw_.value();
        scalar lambda = max(lambda_.value(), SMALL);

        scalar Se = (alphaW - Srw) / (1.0 - Srw);
        Se = max(1e-4, min(1.0 - 1e-4, Se));

        scalar Krw_val, Krnw_val;

        Krw_val = pow(Se, ((2.0 + (3.0 * lambda)) / lambda));
        Krnw_val = pow((1.0 - Se), 2.0) * (1.0 - pow(Se, (2.0 + lambda) / lambda));
  
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

Foam::krBrooksCorey::krBrooksCorey
(
    const fvMesh& mesh,
    const volScalarField& alpha1,
    const volScalarField& alpha2,
    const dictionary& dict
)
:
    relativePermeabilityModelBase(mesh, alpha1, alpha2, dict),
    Srw_("Srw", dimless, dict.getOrDefault<scalar>("Srw", 0.05)),
    lambda_("lambda", dimless, dict.getOrDefault<scalar>("lambda", 2.0))
{
    if (enabled_)
    {
        Info<< "Relative permeability model enabled (Brooks-Corey)" << endl;
        Info<< "    Residual Saturation (Wetting): " << Srw_.value() << endl;
        Info<< "    Lambda: " << lambda_.value() << endl;
        Info<< "    Darcy coefficient: " << D_.value() << " 1/m^2" << endl;
    }
    calculateRelativePermeability();
}

bool Foam::krBrooksCorey::writeData(Ostream& os) const
{
    relativePermeabilityModelBase::writeData(os);
    return os.good();
}

bool Foam::krBrooksCorey::read(const dictionary& dict)
{
    relativePermeabilityModelBase::read(dict);
    Srw_.value() = dict.getOrDefault("Srw", Srw_.value());
    lambda_.value() = dict.getOrDefault("lambda", lambda_.value());
    return true;
}

Foam::Ostream& Foam::operator<<(Foam::Ostream& os, const krBrooksCorey& model)
{
    model.writeData(os);
    return os;
}
