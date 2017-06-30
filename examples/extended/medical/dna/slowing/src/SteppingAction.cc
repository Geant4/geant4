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
// shall cite the following Geant4-DNA collaboration publications:
// Phys. Med. 31 (2015) 861-874
// Med. Phys. 37 (2010) 4692-4708
// The Geant4-DNA web site is available at http://geant4-dna.org
//
// $Id$
//
/// \file medical/dna/slowing/src/SteppingAction.cc
/// \brief Implementation of the SteppingAction class

#include "Analysis.hh"
#include "SteppingAction.hh"
#include "EventAction.hh"

#include "G4SystemOfUnits.hh"
#include "G4Electron.hh"
#include "G4EventManager.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

SteppingAction::SteppingAction(EventAction* event): G4UserSteppingAction()
,fEventAction(event)
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

SteppingAction::~SteppingAction()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void SteppingAction::UserSteppingAction(const G4Step* step)
{ 

 G4double edep = step->GetTotalEnergyDeposit();
 
 //total energy deposit 
 fEventAction->AddEdep(edep);      
 
 if (step->GetTrack()->GetDynamicParticle()->GetDefinition()
     == G4Electron::ElectronDefinition())
 {
   const G4String& processName = step->GetPostStepPoint()->
      GetProcessDefinedStep()->GetProcessName();
   
   if (processName!="Transportation")
   {
      // get analysis manager
      G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();

      // fill histo
    
      // 1) ALL e-
      // analysisManager
      // ->FillH1(1, 
      // step->GetTrack()->GetKineticEnergy()/eV,step->GetStepLength()/nm);
      analysisManager->FillH1(1,
       step->GetPreStepPoint()->GetKineticEnergy()/eV,
       step->GetStepLength()/nm);

      // 2) PARENT e-
      if (step->GetTrack()->GetParentID()==0 
          && step->GetTrack()->GetTrackID()==1) 
      // analysisManager->FillH1(2, step->GetTrack()->GetKineticEnergy()/eV,
      // step->GetStepLength()/nm);
      analysisManager->FillH1(2, 
       step->GetPreStepPoint()->GetKineticEnergy()/eV,
       step->GetStepLength()/nm);
    
      // 3) SECONDARY e-
      if (   !(step->GetTrack()->GetParentID()==0) 
          || !(step->GetTrack()->GetTrackID()==1) ) 
      // analysisManager->FillH1(3, step->GetTrack()->GetKineticEnergy()/eV,
      // step->GetStepLength()/nm);
       analysisManager->FillH1(3, 
        step->GetPreStepPoint()->GetKineticEnergy()/eV,
        step->GetStepLength()/nm);
    
   }
    
 }
}  

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......  
