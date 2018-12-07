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
/// \file electromagnetic/TestEm1/src/TrackingAction.cc
/// \brief Implementation of the TrackingAction class
//
// 
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "TrackingAction.hh"
#include "PrimaryGeneratorAction.hh"
#include "Run.hh"
#include "HistoManager.hh"

#include "G4RunManager.hh"
#include "G4Track.hh"

#include "G4SystemOfUnits.hh"
#include "G4UnitsTable.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

TrackingAction::TrackingAction(PrimaryGeneratorAction* prim)
:G4UserTrackingAction(),fPrimary(prim)
{ }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void TrackingAction::PreUserTrackingAction(const G4Track*)
{
/*
  //debug  
  const G4DynamicParticle* dynamic = aTrack->GetDynamicParticle();
  G4double dynamCharge = dynamic->GetCharge();
  G4int occup          = dynamic->GetTotalOccupancy();
  G4double   pdgMass   = dynamic->GetParticleDefinition()->GetPDGMass();    
  G4double invarMass   = dynamic->Get4Momentum().m();  
  G4double dynamMass   = dynamic->GetMass();
  G4double dif1 = invarMass - pdgMass;
  G4double dif2 = dynamMass - pdgMass;
  
  G4cout
    << "\n  Begin of track :" 
    << "\n    charge= " <<  dynamCharge << "  occupancy= " << occup
    << "\n   pdgMass= " << G4BestUnit (pdgMass  , "Energy")    
///    << "\n invarMass= " << G4BestUnit (invarMass, "Energy")
///    << "   invar-pdg= " << G4BestUnit (dif1     , "Energy")
    << "\n dynamMass= " << G4BestUnit (dynamMass, "Energy")
    << "   dynam-pdg= " << G4BestUnit (dif2     , "Energy")
    << G4endl;          
*/             
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void TrackingAction::PostUserTrackingAction(const G4Track* aTrack)
{
  //increase nb of processed tracks 
  //count nb of steps of this track
  G4int   nbSteps = aTrack->GetCurrentStepNumber();
  G4double Trleng = aTrack->GetTrackLength();
  
  Run* run = static_cast<Run*>(
             G4RunManager::GetRunManager()->GetNonConstCurrentRun()); 
    
  if (aTrack->GetDefinition()->GetPDGCharge() == 0.) {
    run->CountTraks0(1); 
    run->CountSteps0(nbSteps);
  
  } else {
    run->CountTraks1(1); 
    run->CountSteps1(nbSteps);
  }
  
  //true and projected ranges for primary particle
  if (aTrack->GetTrackID() == 1) {
    run->AddTrueRange(Trleng);
    G4ThreeVector vertex = fPrimary->GetParticleGun()->GetParticlePosition();
    G4ThreeVector position = aTrack->GetPosition() - vertex;      
    run->AddProjRange(position.x());
    run->AddTransvDev(position.y());
    run->AddTransvDev(position.z());
    
    G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();
    analysisManager->FillH1(1,Trleng);
    analysisManager->FillH1(2,(float)nbSteps);        
  }
/*
  //debug  
  const G4DynamicParticle* dynamic = aTrack->GetDynamicParticle();
  G4double dynamCharge = dynamic->GetCharge();
  G4int occup          = dynamic->GetTotalOccupancy();
  G4double   pdgMass   = dynamic->GetParticleDefinition()->GetPDGMass();    
  G4double invarMass   = dynamic->Get4Momentum().m();  
  G4double dynamMass   = dynamic->GetMass();
  G4double dif1 = invarMass - pdgMass;
  G4double dif2 = dynamMass - pdgMass;
  
  G4cout
    << "\n  End of track :"    
    << "\n    charge= " <<  dynamCharge << "  occupancy= " << occup
    << "\n   pdgMass= " << G4BestUnit (pdgMass  , "Energy")    
///    << "\n invarMass= " << G4BestUnit (invarMass, "Energy")
///    << "   invar-pdg= " << G4BestUnit (dif1     , "Energy")
    << "\n dynamMass= " << G4BestUnit (dynamMass, "Energy")
    << "   dynam-pdg= " << G4BestUnit (dif2     , "Energy")
    << G4endl;          
*/                       
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

