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
/// \file electromagnetic/TestEm3/src/EventAction.cc
/// \brief Implementation of the EventAction class
//
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "EventAction.hh"

#include "Run.hh"
#include "HistoManager.hh"

#include "G4RunManager.hh"
#include "G4Event.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

EventAction::EventAction(DetectorConstruction* det)
:fDetector(det)
{ }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void EventAction::BeginOfEventAction(const G4Event*)
{       
  //initialize EnergyDeposit per event
  //
  for (G4int k=0; k<kMaxAbsor; k++) {
    fEnergyDeposit[k] = fTrackLengthCh[k] = 0.0;   
  }
  // initialize EnergyLeakage per event
  //
  fEnergyLeak = 0.0;    
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void EventAction::SumEnergy(G4int k, G4double de, G4double dl)
{       
  fEnergyDeposit[k] += de;
  fTrackLengthCh[k] += dl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void EventAction::SumEnergyLeak(G4double eleak)
{       
  fEnergyLeak += eleak;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void EventAction::EndOfEventAction(const G4Event*)
{
  //get Run
  Run* run = static_cast<Run*>(
             G4RunManager::GetRunManager()->GetNonConstCurrentRun());
 
  G4double EdepTot =0.;            
  for (G4int k=1; k<=fDetector->GetNbOfAbsor(); k++) {
     run->FillPerEvent(k,fEnergyDeposit[k],fTrackLengthCh[k]);
     if (fEnergyDeposit[k] > 0.)
             G4AnalysisManager::Instance()->FillH1(k, fEnergyDeposit[k]);
     EdepTot += fEnergyDeposit[k];	     
  }
  
  run->SumEnergies(EdepTot,fEnergyLeak);
  
  //histograms
  G4AnalysisManager* analysis = G4AnalysisManager::Instance();
  analysis->FillH1(kMaxAbsor, EdepTot);
  G4int id = 2*kMaxAbsor+3;
  analysis->FillH1(id, fEnergyLeak);
  G4double ETot = EdepTot + fEnergyLeak;
  analysis->FillH1(++id, ETot); 	 
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

