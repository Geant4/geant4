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
// SteppingAction.cc
// 
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "SteppingAction.hh"

#include "DetectorConstruction.hh"
#include "RunAction.hh"
#include "EventAction.hh"
#include "StackingAction.hh"

#include "G4ProcessManager.hh"
#include "G4Step.hh"

using namespace std; 
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

SteppingAction::SteppingAction(DetectorConstruction* DET, RunAction* RA,
                               EventAction* EA )
:detector(DET), runaction(RA), eventaction(EA) 
{ }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

SteppingAction::~SteppingAction()
{ }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void SteppingAction::UserSteppingAction(const G4Step* aStep )
{

	//energy loss in detector
 if (aStep->GetPreStepPoint()->GetTouchableHandle()->GetVolume() 
     != detector->GetAbsorber() ||
     aStep->GetPostStepPoint()->GetTouchableHandle()->GetVolume()
     != detector->GetAbsorber()) return;
 

  
 G4int IDp =  aStep->GetTrack()->GetParentID();

	//primary particle
 if (IDp==0 ){


 	eventaction->AddEnergy(aStep->GetTotalEnergyDeposit());
 	eventaction->AddNonIonEnergy (aStep->GetNonIonizingEnergyDeposit());


 	eventaction->AddTrakLenPrim(aStep->GetStepLength());
 	eventaction->CountStepsPrim();

 	} 


}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

