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
#include "ActionInitialization.hh"
#include "G4SystemOfUnits.hh"
#include "SteppingAction.hh"
#include "RunAction.hh"
#include "StackingAction.hh"
#include "DetectorConstruction.hh"
#include "G4Scheduler.hh"
#include "TimeStepAction.hh"
#include "ITTrackingInteractivity.hh"
#include "G4RunManager.hh"
#include "G4MoleculeGun.hh"
#include "PrimaryGeneratorAction.hh"
#include "G4DNAChemistryManager.hh"
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

ActionInitialization::ActionInitialization(DetectorConstruction* pDetector)
    : G4VUserActionInitialization()
    , fpDetector(pDetector)
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

ActionInitialization::~ActionInitialization()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void ActionInitialization::BuildForMaster() const
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void ActionInitialization::Build() const
{
    SetUserAction(new PrimaryGeneratorAction);
    RunAction* pRunAction = new RunAction();
    SetUserAction(pRunAction);
    SteppingAction* pSteppingAction = new SteppingAction(fpDetector);
    SetUserAction(pSteppingAction);
    SetUserAction(new StackingAction());
    if(G4DNAChemistryManager::IsActivated())
    {
        G4Scheduler::Instance()->
        SetUserAction(new TimeStepAction());
//stop at this time
        G4Scheduler::Instance()->SetEndTime(2.5*nanosecond);
        G4Scheduler::Instance()->SetVerbose(1);
        ITTrackingInteractivity* itInteractivity = 
        new ITTrackingInteractivity();
        G4Scheduler::Instance()->SetInteractivity(itInteractivity);
        G4DNAChemistryManager::Instance()->
        SetGun(((DetectorConstruction*)fpDetector)->GetGun());
    }
}
