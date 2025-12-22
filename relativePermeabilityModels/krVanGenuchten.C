/*---------------------------------------------------------------------------*\
  =========                 |\
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |\
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |\
-------------------------------------------------------------------------------
    Copyright (C) 2024 Custom Implementation
\*---------------------------------------------------------------------------*/

#include "krVanGenuchten.H"
#include "addToRunTimeSelectionTable.H"
#include "fvc.H"
#include "IOstream.H"

namespace Foam
{
    defineTypeNameAndDebug(krVanGenuchten, 0);
    addToRunTimeSelectionTable(relativePermeabilityModelBase, krVanGenuchten, dictionary);
}

void Foam::krVanGenuchten::calculateRelativePermeability()
{
    forAll(Kr1_, cellI)
    {
        scalar alphaW = (wettingPhase_ == 1) ? alpha1_[cellI] : alpha2_[cellI];

        scalar Srw = Srw_.value();
        scalar m = max(m_.value(), SMALL);

        scalar Se = (alphaW - Srw) / (1.0 - Srw);
        Se = max(1e-4, min(1.0 - 1e-4, Se));

        scalar Krw_val, Krnw_val;

        Krw_val = pow(Se, 0.5) * pow(1.0 - pow((1.0 - pow(Se, (1.0/m))), m), 2.0);
        Krnw_val = pow((1.0 - Se), 0.5) * pow((1.0 - pow(Se, (1.0/m))), (2.0*m));

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

Foam::krVanGenuchten::krVanGenuchten
(
    const fvMesh& mesh,
    const volScalarField& alpha1,
    const volScalarField& alpha2,
    const dictionary& dict
)
:
    relativePermeabilityModelBase(mesh, alpha1, alpha2, dict),
    Srw_("Srw", dimless, dict.getOrDefault<scalar>("Srw", 0.05)),
    m_("m", dimless, dict.getOrDefault<scalar>("m", 0.5))
{
    if (enabled_)
    {
        Info<< "Relative permeability model enabled (Van Genuchten)" << endl;
        Info<< "    Residual Saturation (Wetting): " << Srw_.value() << endl;
        Info<< "    m: " << m_.value() << endl;
        Info<< "    Darcy coefficient: " << D_.value() << " 1/m^2" << endl;
    }
    calculateRelativePermeability();
}

bool Foam::krVanGenuchten::writeData(Ostream& os) const
{
    relativePermeabilityModelBase::writeData(os);
    return os.good();
}

bool Foam::krVanGenuchten::read(const dictionary& dict)
{
    relativePermeabilityModelBase::read(dict);
    Srw_.value() = dict.getOrDefault("Srw", Srw_.value());
    m_.value() = dict.getOrDefault("m", m_.value());
    return true;
}

Foam::Ostream& Foam::operator<<(Foam::Ostream& os, const krVanGenuchten& model)
{
    model.writeData(os);
    return os;
}
