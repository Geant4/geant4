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
/// \file CalorimeterSD.cc
/// \brief Implementation of the CalorimeterSD class

#include "Analysis.hh" 

#include "CalorimeterSD.hh"
#include "DetectorParameters.hh"
#include "ApplicationParameters.hh"
#include "G4HCofThisEvent.hh"
#include "G4Step.hh"
#include "G4ThreeVector.hh"
#include "G4SDManager.hh"
#include "G4SystemOfUnits.hh"
#include "G4ios.hh"

#include <cmath>

using namespace ApplicationParameters;
using namespace DetectorParameters;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

CalorimeterSD::CalorimeterSD(const G4String& name)
 : G4VSensitiveDetector(name),
   fZHitsCollection(0),
   fXYHitsCollection(0)
{
  collectionName.insert("zHitsCollection");
  collectionName.insert("xyHitsCollection");
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

CalorimeterSD::~CalorimeterSD() 
{ 
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void CalorimeterSD::Initialize(G4HCofThisEvent* hce)
{
  // Create hits collection
  fZHitsCollection 
    = new CalorHitsCollection(SensitiveDetectorName, collectionName[0]); 
  fXYHitsCollection 
    = new CalorHitsCollection(SensitiveDetectorName, collectionName[1]); 

  // Add this collection in hce
  G4int zhcID 
    = G4SDManager::GetSDMpointer()->GetCollectionID(collectionName[0]);
  hce->AddHitsCollection( zhcID, fZHitsCollection ); 
  
  G4int xyhcID 
    = G4SDManager::GetSDMpointer()->GetCollectionID(collectionName[1]);
  hce->AddHitsCollection( xyhcID, fXYHitsCollection ); 
  
  // Create hits
  // fNofCells for cells + one more for total sums 
  for (G4int i=0; i<NofZLayers+1; i++ ) {
    fZHitsCollection->insert(new CalorHit());
  }
  for (G4int i=0; i<NofXYLayers*NofXYLayers+1; i++ ) {
    fXYHitsCollection->insert(new CalorHit());
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4bool CalorimeterSD::ProcessHits(G4Step* step, 
                                     G4TouchableHistory*)
{  
  // energy deposit
  G4double edep = step->GetTotalEnergyDeposit();
  
  // step length
  G4double stepLength = 0.;
  if ( step->GetTrack()->GetDefinition()->GetPDGCharge() != 0. ) {
    stepLength = step->GetStepLength();
  }

  if ( edep==0. && stepLength == 0. ) return false;      

  G4TouchableHistory* touchable
    = (G4TouchableHistory*)(step->GetPreStepPoint()->GetTouchable());
    
  // Get calorimeter cell id 
  // G4cout << "got touchable " <<  touchable->GetVolume(2)->GetName() << G4endl;
  G4int zLayerNumber = touchable->GetReplicaNumber(2);
  G4int yLayerNumber = touchable->GetReplicaNumber(1);
  G4int xLayerNumber = touchable->GetReplicaNumber();
  G4int xyLayerNumber 
    = xLayerNumber*NofXYLayers + yLayerNumber;
  
  // Get hit accounting data for this cell
  // Get hit accounting data for this cell
  CalorHit* zhit = (*fZHitsCollection)[zLayerNumber];
  if ( ! zhit ) {
    G4ExceptionDescription msg;
    msg << "Cannot access hit " << zLayerNumber; 
    G4Exception("CalorimeterSD::ProcessHits()",
      "MyCode0004", FatalException, msg);
  }         

  // Get hit accounting data for this cell
  CalorHit* xyhit = (*fXYHitsCollection)[xyLayerNumber];
  if ( ! xyhit ) {
    G4ExceptionDescription msg;
    msg << "Cannot access hit " << xyLayerNumber; 
    G4Exception("CalorimeterSD::ProcessHits()",
      "MyCode0004", FatalException, msg);
  }         

  // Get hit for total accounting
  CalorHit* zhitTotal 
    = (*fZHitsCollection)[fZHitsCollection->entries()-1];
  CalorHit* xyhitTotal 
    = (*fXYHitsCollection)[fXYHitsCollection->entries()-1];
  
  // Add values
  zhit->Add(edep, stepLength);
  xyhit->Add(edep, stepLength);
  zhitTotal->Add(edep, stepLength); 
  xyhitTotal->Add(edep, stepLength);

  // Fill histograms
  G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();

  // Momentum in z layer #2
  if ( zLayerNumber == 2 ) {
    G4ThreeVector momentum = step->GetTrack()->GetMomentum();
    if ( TestH1 ) analysisManager->FillH1(10, momentum.x());
    if ( TestH1 ) analysisManager->FillH1(11, momentum.x());
    if ( TestH2 ) analysisManager->FillH2(0, momentum.x(), momentum.y());
    if ( TestH2 ) analysisManager->FillH2(1, momentum.x(), momentum.y());
    if ( TestH3 ) analysisManager->FillH3(0, momentum.x(), momentum.y(), momentum.z());
    if ( TestH3 ) analysisManager->FillH3(1, momentum.x(), momentum.y(), momentum.z());
  }  
         
  return true;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void CalorimeterSD::EndOfEvent(G4HCofThisEvent*)
{
  G4double dLCumul = 0.;
  G4double sumEdep = 0.;  

  // Fill 1D profiles
  G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();
  for ( G4int i=0; i<NofZLayers; i++ ) {
    sumEdep +=  (*fZHitsCollection)[i]->GetEdep(); 
    dLCumul = (i + 0.5) * LayerZSize;
    if ( TestP1 ) analysisManager->FillP1(0, dLCumul, sumEdep/dLCumul);
    if ( TestP1 ) analysisManager->FillP1(1, dLCumul, sumEdep/dLCumul);
  }  

  G4double dXCumul = 0.;
  G4double dYCumul = 0.;

  // Fill 2D profiles
  for ( G4int ix=0; ix< NofXYLayers; ix++ ) {
    dXCumul = (ix + 0.5 - NofXYLayers*0.5) * LayerXYSize;
    for ( G4int iy=0; iy< NofXYLayers; iy++ ) {
      dYCumul = (iy + 0.5 - NofXYLayers*0.5) * LayerXYSize;
      G4int ihit = ix*NofXYLayers + iy;
      G4double edep = (*fXYHitsCollection)[ihit]->GetEdep(); 
      if ( edep == 0 ) continue;
      if ( TestP2 ) analysisManager->FillP2(0, dXCumul, dYCumul, edep);
      if ( TestP2 ) analysisManager->FillP2(1, dXCumul, dYCumul, edep);
    }    
  }  

  G4int nofZHits = fZHitsCollection->entries();
  if ( verboseLevel>1 ) { 
     G4cout << "\n-------->Z Hits Collection: in this event they are " << nofZHits 
            << " hits: " << G4endl;
     for ( G4int i=0; i<nofZHits; i++ ) (*fZHitsCollection)[i]->Print();
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
