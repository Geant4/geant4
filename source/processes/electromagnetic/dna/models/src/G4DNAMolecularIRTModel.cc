//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
//
/*
 * G4DNAMolecularIRTModel.cc
 *
 *  Created on: Jul 23, 2019
 *      Author: W. G. Shin
 *              J. Ramos-Mendez and B. Faddegon
*/

#include <globals.hh>
#include <G4DNAMolecularReactionTable.hh>
#include <G4DNAMolecularIRTModel.hh>
#include <G4DNASmoluchowskiReactionModel.hh>
#include <G4ExceptionSeverity.hh>
#include <G4Molecule.hh>
#include <G4ReferenceCast.hh>

#include "G4DNAIRT.hh"
#include "G4DNAIRTMoleculeEncounterStepper.hh"

G4DNAMolecularIRTModel::G4DNAMolecularIRTModel(const G4String& name)
    : G4DNAMolecularIRTModel(name,
                std::unique_ptr<G4DNAIRTMoleculeEncounterStepper>(new G4DNAIRTMoleculeEncounterStepper()),
                std::unique_ptr<G4DNAIRT>(new G4DNAIRT()))
{
}

G4DNAMolecularIRTModel::G4DNAMolecularIRTModel(const G4String& name,
                                                             std::unique_ptr<G4VITTimeStepComputer> pTimeStepper,
                                                             std::unique_ptr<G4VITReactionProcess> pReactionProcess)
    : G4VITStepModel(std::move(pTimeStepper),
                     std::move(pReactionProcess),
                     name)
    , fMolecularReactionTable(reference_cast<const G4DNAMolecularReactionTable*>(fpReactionTable))
{
    fType1 = G4Molecule::ITType();
    fType2 = G4Molecule::ITType();
}

G4DNAMolecularIRTModel::~G4DNAMolecularIRTModel() = default;

void G4DNAMolecularIRTModel::Initialize()
{
    if(fpReactionTable == nullptr)
    {
        SetReactionTable(G4DNAMolecularReactionTable::GetReactionTable());
    }

    if(!fpReactionModel)
    {
        fpReactionModel.reset(new G4DNASmoluchowskiReactionModel());
    }

    fpReactionModel->SetReactionTable((const G4DNAMolecularReactionTable*) fpReactionTable);

    ((G4DNAIRT*) fpReactionProcess.get())->SetReactionModel(fpReactionModel.get());
    ((G4DNAIRTMoleculeEncounterStepper*) fpTimeStepper.get())->SetReactionModel(fpReactionModel.get());

    G4VITStepModel::Initialize();
}

void G4DNAMolecularIRTModel::PrintInfo()
{
#ifdef G4VERBOSE
    G4cout << fName << " will be used" << G4endl;
#endif
}

void G4DNAMolecularIRTModel::SetReactionModel(G4VDNAReactionModel* pReactionModel)
{
    fpReactionModel.reset(pReactionModel);
}

G4VDNAReactionModel* G4DNAMolecularIRTModel::GetReactionModel()
{
    return fpReactionModel.get();
}
