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
// This example is provided by the Geant4-DNA collaboration
// Any report or published results obtained using the Geant4-DNA software
// shall cite the following Geant4-DNA collaboration publication:
// Med. Phys. 37 (2010) 4692-4708
// J. Comput. Phys. 274 (2014) 841-882
// The Geant4-DNA web site is available at http://geant4-dna.org
//
//
/// \file DetectorConstruction.cc
/// \brief Implementation of the DetectorConstruction class

#include "ActionInitialization.hh"

#include "PrimaryGeneratorAction.hh"
#include "RunAction.hh"
#include "SteppingAction.hh"
#include "DetectorConstruction.hh"
#include "TrackingAction.hh"
#include "G4RunManager.hh"

#include "PrimaryGeneratorAction.hh"
#include "RunAction.hh"
#include "SteppingAction.hh"

#include "G4DNAChemistryManager.hh"

#include "G4SystemOfUnits.hh"
#include "G4UnitsTable.hh"

#include "G4Threading.hh"

#include "TrackingAction.hh"

// chemistry
#include "G4Scheduler.hh"
#include "StackingAction.hh"
#include "TimeStepAction.hh"
#include "ITTrackingInteractivity.hh"
#include "ITSteppingAction.hh"
#include "ITTrackingAction.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

ActionInitialization::ActionInitialization() : G4VUserActionInitialization()
{
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

ActionInitialization::~ActionInitialization()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void ActionInitialization::BuildForMaster() const
{
 // In MT mode, to be clearer, the RunAction class for the master thread might
 // be different than the one used for the workers.
 // This RunAction will be called before and after starting the
 // workers.
 // For more details, please refer to :
 // https://twiki.cern.ch/twiki/bin/view/Geant4/Geant4MTForApplicationDevelopers
 //
 // RunAction* runAction= new RunAction();
 // SetUserAction(runAction);
}

void ActionInitialization::Build() const
{
  PrimaryGeneratorAction* primGenAction = new PrimaryGeneratorAction;
  SetUserAction(primGenAction);

  // Set optional user action classes
  SetUserAction(new RunAction());
  SetUserAction(new TrackingAction());
  SetUserAction(new SteppingAction());
  SetUserAction(new StackingAction());

  // chemistry part
  if(G4DNAChemistryManager::IsActivated()){
    G4Scheduler::Instance()->SetUserAction(new TimeStepAction());

    // Uncomment and set to stop chemistry stage after:
    // ...given number of time steps
    //G4Scheduler::Instance()->SetMaxNbSteps(1000);

    // ...OR reaching this time
    G4Scheduler::Instance()->SetEndTime(100*nanosecond);

    G4Scheduler::Instance()->SetVerbose(1);

    ITTrackingInteractivity* itInteractivity = new ITTrackingInteractivity();
    itInteractivity->SetUserAction(new ITSteppingAction);
    itInteractivity->SetUserAction(new ITTrackingAction);
    G4Scheduler::Instance()->SetInteractivity(itInteractivity);
  }
/*
  // To output the pre-chemical stage
  //
  G4String fileName ("output");

  if(G4RunManager::GetRunManager()->GetRunManagerType() ==
      G4RunManager::sequentialRM)
  {
    // write initial situation at 1 picosecond
    G4DNAChemistryManager::Instance()->WriteInto(fileName + ".txt");
  }
  else
  {
    G4int id = G4Threading::G4GetThreadId();

    G4String fileName_mt = fileName;
    fileName_mt += G4UIcommand::ConvertToString(id);
    fileName_mt += ".txt";

    G4cout << "chosen file name : " << fileName_mt << G4endl;

    G4DNAChemistryManager::Instance()->WriteInto(fileName_mt);
  }
*/
}  
