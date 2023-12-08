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
/// \file ChemActionInitialization.cc
/// \brief Implementation of the ChemActionInitialization class


#include "ChemActionInitialization.hh"
#include "ChemNtupleManager.hh"
#include "ChemRunAction.hh"

#include "G4Timer.hh"
#include "G4UnitsTable.hh"

#include "ChemStackingAction.hh"
#include "ChemPrimaryGeneratorAction.hh"

// chemistry
#include "G4Scheduler.hh"
#include "G4DNAChemistryManager.hh"
#include "ChemITSteppingAction.hh"
#include "ChemTimeStepAction.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

ChemActionInitialization::ChemActionInitialization(ChemNtupleManager* nMana, ChemPhysicsList *phys)
:fpNtuple(nMana), fPhysList(phys)
{

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void ChemActionInitialization::BuildForMaster() const
{
    ;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void ChemActionInitialization::Build() const
{
    SetUserAction(new ChemRunAction(fpNtuple));
    
    SetUserAction(new ChemPrimaryGeneratorAction());

    SetUserAction(new ChemStackingAction());

    G4bool chemistryFlag = G4DNAChemistryManager::Instance()->IsActivated();

    if(chemistryFlag)
    {
        G4Scheduler::Instance()->SetVerbose(0);

        G4Scheduler::Instance()->SetMaxZeroTimeAllowed(10000);
    
        ChemTimeStepAction* timeStepAction = new ChemTimeStepAction(fpNtuple,fPhysList->GetTimeStepModel());
        G4Scheduler::Instance()->SetUserAction(timeStepAction); 
    }
    
}  

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......