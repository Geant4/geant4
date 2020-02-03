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

#include <globals.hh>
#include <G4DNAMolecularReaction.hh>
#include <G4DNAMolecularReactionTable.hh>
#include <G4DNAMolecularStepByStepModel.hh>
#include <G4DNAMoleculeEncounterStepper.hh>
#include <G4DNASmoluchowskiReactionModel.hh>
#include <G4ExceptionSeverity.hh>
#include <G4Molecule.hh>
#include <G4ReferenceCast.hh>
#include <G4VDNAReactionModel.hh>

G4DNAMolecularStepByStepModel::G4DNAMolecularStepByStepModel(const G4String& name)
    : G4DNAMolecularStepByStepModel(name,
                                    std::unique_ptr<G4DNAMoleculeEncounterStepper>(new G4DNAMoleculeEncounterStepper()),
                                    std::unique_ptr<G4DNAMolecularReaction>(new G4DNAMolecularReaction()))
{
}

G4DNAMolecularStepByStepModel::G4DNAMolecularStepByStepModel(const G4String& name,
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

G4DNAMolecularStepByStepModel::~G4DNAMolecularStepByStepModel() = default;

void G4DNAMolecularStepByStepModel::Initialize()
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

    ((G4DNAMolecularReaction*) fpReactionProcess.get())->SetReactionModel(fpReactionModel.get());
    ((G4DNAMoleculeEncounterStepper*) fpTimeStepper.get())->SetReactionModel(fpReactionModel.get());

    G4VITStepModel::Initialize();
}

void G4DNAMolecularStepByStepModel::PrintInfo()
{
#ifdef G4VERBOSE
    G4cout << fName << " will be used" << G4endl;
#endif
}

void G4DNAMolecularStepByStepModel::SetReactionModel(G4VDNAReactionModel* pReactionModel)
{
    fpReactionModel.reset(pReactionModel);
}

G4VDNAReactionModel* G4DNAMolecularStepByStepModel::GetReactionModel()
{
    return fpReactionModel.get();
}
