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
// $Id: G4DNAMolecularStepByStepModel.cc 65022 2012-11-12 16:43:12Z gcosmo $
//
#include "G4DNAMolecularStepByStepModel.hh"
#include "G4VDNAReactionModel.hh"

G4DNAMolecularStepByStepModel::G4DNAMolecularStepByStepModel(const G4String& name) :
    G4VITModel(name),
    fMolecularReactionTable(reference_cast<const G4DNAMolecularReactionTable*>(fpReactionTable))
{
    fpTimeStepper = new G4DNAMoleculeEncounterStepper();
    fpReactionProcess = new G4DNAMolecularReaction();

    fType1 = G4Molecule::ITType();
    fType2 = G4Molecule::ITType();
    fReactionModel = 0;
}

G4DNAMolecularStepByStepModel::~G4DNAMolecularStepByStepModel()
{
    if(fReactionModel) delete fReactionModel;
}

G4DNAMolecularStepByStepModel& G4DNAMolecularStepByStepModel::operator=(const G4DNAMolecularStepByStepModel& right)
{
    G4ExceptionDescription exceptionDescription
            ("Use copy constructor rather than assignement operator.");
    G4Exception("G4DNAMolecularStepByStepModel::operator=(const G4DNAMolecularStepByStepModel&)","G4DNAMolecularStepByStepModel001",
                FatalErrorInArgument,exceptionDescription);

    if (&right==this) return *this;
    return *this; // avoid warnings
}


G4DNAMolecularStepByStepModel::G4DNAMolecularStepByStepModel(const G4DNAMolecularStepByStepModel& right):
    G4VITModel(right),
    fMolecularReactionTable(reference_cast<const G4DNAMolecularReactionTable*>(fpReactionTable))
{
        fpReactionTable = right.fpReactionTable;
        if(right.fReactionModel)
        {
            fReactionModel = right.fReactionModel->Clone();
            ((G4DNAMolecularReaction*)  fpReactionProcess)->SetReactionModel(fReactionModel);
            ((G4DNAMoleculeEncounterStepper*) 	fpTimeStepper)->SetReactionModel(fReactionModel);
        }
        else fReactionModel = 0;
}

void G4DNAMolecularStepByStepModel::Initialize()
{
    fReactionModel ->SetReactionTable((const G4DNAMolecularReactionTable*)fpReactionTable);
    G4VITModel::Initialize();
}

void G4DNAMolecularStepByStepModel::PrintInfo()
{
#ifdef G4VERBOSE
    G4cout << "DNAMolecularStepByStepModel will be used" << G4endl;
#endif
}
