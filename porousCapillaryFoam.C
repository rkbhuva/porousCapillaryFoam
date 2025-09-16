/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
    Copyright (C) 2020 OpenCFD Ltd.
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM. If not, see <http://www.gnu.org/licenses/>.

Application
    twoPhaseEulerFoam

Group
    grpMultiphaseSolvers

Description
    Solver for a system of two compressible fluid phases with one dispersed
    phase. Eg, gas bubbles in a liquid including heat-transfer.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "twoPhaseSystem.H"
#include "PhaseCompressibleTurbulenceModel.H"
#include "pimpleControl.H"
#include "fvOptions.H"
#include "fixedValueFvsPatchFields.H"
#include "capillaryPressureModelBase.H"
#include "relativePermeabilityModelBase.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    argList::addNote
    (
        "Solver for porous two phase flow isothermal. \n"
    );

    #include "postProcess.H"

    #include "addCheckCaseOptions.H"
    #include "setRootCaseLists.H"
    #include "createTime.H"
    #include "createMesh.H"
    #include "createControl.H"
    #include "createFields.H"
    #include "createFieldRefs.H"

    Info<< "Reading capillaryTransportProperties\n" << endl;
    IOdictionary capillaryTransportModelDict
    (
        IOobject
        (
            "capillaryTransportProperties",
            runTime.constant(),
            mesh,
            IOobject::MUST_READ_IF_MODIFIED,
            IOobject::NO_WRITE
        )
    );

    Info<< "Creating capillaryPressureModel and relativePermeabilityModel\n" << endl;
    
    autoPtr<Foam::capillaryPressureModelBase> capPressureModel
    (
        Foam::capillaryPressureModelBase::New
        (
            mesh, alpha1, alpha2, capillaryTransportModelDict
        )
    );

    autoPtr<Foam::relativePermeabilityModelBase> relPermModel
    (
        Foam::relativePermeabilityModelBase::New
        (
            mesh, alpha1, alpha2, capillaryTransportModelDict
        )
    );

    #include "createTimeControls.H"
    #include "CourantNos.H"
    #include "setInitialDeltaT.H"

    bool faceMomentum
    (
        pimple.dict().getOrDefault("faceMomentum", false)
    );

    bool implicitPhasePressure
    (
        mesh.solverDict(alpha1.name()).getOrDefault
        (
            "implicitPhasePressure", false
        )
    );

    #include "pUf/createDDtU.H"
    #include "pU/createDDtU.H"

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< "\nStarting time loop\n" << endl;

    while (runTime.run())
    {
        #include "readTimeControls.H"
        #include "CourantNos.H"
        #include "setDeltaT.H"

        ++runTime;
        Info<< "Time = " << runTime.timeName() << nl << endl;

        // --- Pressure-velocity PIMPLE corrector loop
        while (pimple.loop())
        {
            fluid.solve();
            fluid.correct();

            #include "contErrs.H"

            if (!capPressureModel.valid())
            {
                FatalErrorInFunction
                    << "capPressureModel is not valid. Check 'capillaryPressureModel' in constant/capillaryTransportProperties."
                    << abort(FatalError);
            }
            if (!relPermModel.valid())
            {
                FatalErrorInFunction
                    << "relPermModel is not valid. Check 'relativePermeabilityModel' in constant/capillaryTransportProperties."
                    << abort(FatalError);
            }

            capPressureModel->update();
            relPermModel->update();

            // Diagnostics
            Info << "Capillary Pressure - min: "
                 << min(capPressureModel->pCapillary()).value()
                 << ", max: " << max(capPressureModel->pCapillary()).value()
                 << ", mean: " << capPressureModel->pCapillary().average().value()
                 << endl;

            Info << "Relative Permeability Phase1 (Kr1) - min: "
                 << min(relPermModel->Kr1()).value()
                 << ", max: " << max(relPermModel->Kr1()).value()
                 << ", mean: " << relPermModel->Kr1().average().value()
                 << endl;

            Info << "Relative Permeability Phase2 (Kr2) - min: "
                 << min(relPermModel->Kr2()).value()
                 << ", max: " << max(relPermModel->Kr2()).value()
                 << ", mean: " << relPermModel->Kr2().average().value()
                 << endl;

            if (faceMomentum)
            {
                #include "pUf/UEqns.H"
                //#include "EEqns.H"
                #include "pUf/pEqn.H"
                #include "pUf/DDtU.H"
            }
            else
            {
                #include "pU/UEqns.H"
                //#include "EEqns.H"
                #include "pU/pEqn.H"
                #include "pU/DDtU.H"
            }

            if (pimple.turbCorr())
            {
                fluid.correctTurbulence();
            }
        }

        #include "write.H"

        Info << "min/max U1: " << min(U1).value() << ", "<< max(U1).value() << endl;
        Info << "min/max U2: " << min(U2).value() << ", "<< max(U2).value() << endl;
        Info << "average alpha1/alpha2: " << alpha1.average().value() << ", "
        << alpha2.average().value() << endl;

        runTime.printExecutionTime(Info);
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
