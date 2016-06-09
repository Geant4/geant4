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
// $Id$
// 
/// \file B4cEventAction.cc
/// \brief Implementation of the B4cEventAction class

#include "B4cEventAction.hh"
#include "B4cCalorimeterSD.hh"
#include "B4cCalorHit.hh"
#include "B4Analysis.hh"

#include "G4RunManager.hh"
#include "G4Event.hh"
#include "G4SDManager.hh"
#include "G4HCofThisEvent.hh"
#include "G4GenericMessenger.hh"
#include "G4UnitsTable.hh"

#include "Randomize.hh"
#include <iomanip>

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

B4cEventAction::B4cEventAction()
 : G4UserEventAction(),
   fMessenger(0),
   fPrintModulo(1)
{
  // Define /B4/event commands using generic messenger class
  fMessenger = new G4GenericMessenger(this, "/B4/event/", "Event control");

  // Define /B4/event/setPrintModulo command
  G4GenericMessenger::Command& setPrintModulo 
    = fMessenger->DeclareProperty("setPrintModulo", 
                                  fPrintModulo, 
                                 "Print events modulo n");
  setPrintModulo.SetRange("value>0");                                
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

B4cEventAction::~B4cEventAction()
{
  delete fMessenger;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

B4cCalorHitsCollection* 
B4cEventAction::GetHitsCollection(const G4String& hcName,
                                  const G4Event* event) const
{
  G4int hcID 
    = G4SDManager::GetSDMpointer()->GetCollectionID(hcName);
  B4cCalorHitsCollection* hitsCollection 
    = static_cast<B4cCalorHitsCollection*>(
        event->GetHCofThisEvent()->GetHC(hcID));
  
  if ( ! hitsCollection ) {
    G4cerr << "Cannot access hitsCollection " << hcName << G4endl;
    exit(1);
  }         

  return hitsCollection;
}    

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void B4cEventAction::PrintEventStatistics(
                              G4double absoEdep, G4double absoTrackLength,
                              G4double gapEdep, G4double gapTrackLength) const
{
  // print event statistics
  G4cout
     << "   Absorber: total energy: " 
     << std::setw(7) << G4BestUnit(absoEdep, "Energy")
     << "       total track length: " 
     << std::setw(7) << G4BestUnit(absoTrackLength, "Length")
     << G4endl
     << "        Gap: total energy: " 
     << std::setw(7) << G4BestUnit(gapEdep, "Energy")
     << "       total track length: " 
     << std::setw(7) << G4BestUnit(gapTrackLength, "Length")
     << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void B4cEventAction::BeginOfEventAction(const G4Event* event)
{  

  G4int eventID = event->GetEventID();
  if ( eventID % fPrintModulo == 0 )  { 
    G4cout << "\n---> Begin of event: " << eventID << G4endl;
    //CLHEP::HepRandom::showEngineStatus();
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void B4cEventAction::EndOfEventAction(const G4Event* event)
{  
  // Get hits collections
  B4cCalorHitsCollection* absoHC
    = GetHitsCollection("AbsorberHitsCollection", event);
  B4cCalorHitsCollection* gapHC
    = GetHitsCollection("GapHitsCollection", event);

  // Get hit with total values
  B4cCalorHit* absoHit = (*absoHC)[absoHC->entries()-1];
  B4cCalorHit* gapHit = (*gapHC)[absoHC->entries()-1];
 
  // Print per event (modulo n)
  //
  G4int eventID = event->GetEventID();
  if ( eventID % fPrintModulo == 0) {
    G4cout << "---> End of event: " << eventID << G4endl;     

    PrintEventStatistics(
      absoHit->GetEdep(), absoHit->GetTrackLength(),
      gapHit->GetEdep(), gapHit->GetTrackLength());
  }  
  
  // Fill histograms, ntuple
  //

  // get analysis manager
  G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();
 
  // fill histograms
  analysisManager->FillH1(1, absoHit->GetEdep());
  analysisManager->FillH1(2, gapHit->GetEdep());
  analysisManager->FillH1(3, absoHit->GetTrackLength());
  analysisManager->FillH1(4, gapHit->GetTrackLength());
  
  // fill ntuple
  analysisManager->FillNtupleDColumn(0, absoHit->GetEdep());
  analysisManager->FillNtupleDColumn(1, gapHit->GetEdep());
  analysisManager->FillNtupleDColumn(2, absoHit->GetTrackLength());
  analysisManager->FillNtupleDColumn(3, gapHit->GetTrackLength());
  analysisManager->AddNtupleRow();  
}  

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
