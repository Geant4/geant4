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
// Authors: Susanna Guatelli and Francesco Romano
// susanna@uow.edu.au, francesco.romano@ct.infn.it
//
//  code based on the basic example B2
//
#include "SensitiveDetector.hh"
#include "G4HCofThisEvent.hh"
#include "G4Step.hh"
#include "G4ThreeVector.hh"
#include "G4SDManager.hh"
#include "G4ios.hh"
#include "G4SystemOfUnits.hh"
#include "G4RunManager.hh"

SensitiveDetector::SensitiveDetector(const G4String& name,
                         const G4String& hitsCollectionName, 
                         G4bool isThisAMicrodosimeter,
                         AnalysisManager* analysis_manager) 
 : G4VSensitiveDetector(name),
   fHitsCollection(NULL)
{
  collectionName.insert(hitsCollectionName);
  analysis = analysis_manager;
  isMicrod = isThisAMicrodosimeter; // only false for the E-stage of the telescope
}

SensitiveDetector::~SensitiveDetector() 
{}

void SensitiveDetector::Initialize(G4HCofThisEvent* hce)
{
  // Create hits collection

  fHitsCollection 
    = new SensitiveDetectorHitsCollection(SensitiveDetectorName, collectionName[0]); 

  // Add this collection in hce

  G4int hcID 
    = G4SDManager::GetSDMpointer()->GetCollectionID(collectionName[0]);
  hce->AddHitsCollection( hcID, fHitsCollection ); 
}

G4bool SensitiveDetector::ProcessHits(G4Step* aStep, 
                                     G4TouchableHistory*)
{  
  // energy deposit
  G4double edep = aStep->GetTotalEnergyDeposit();
  
  if (edep==0.) return false;

  // Uncomment the following lines to only score specific physical volumes 
  /*
  G4String volumeName = aStep -> GetPreStepPoint() -> GetPhysicalVolume()-> GetName();
  if(volumeName != "SV_phys1" && volumeName != "sen_bridge" && volumeName != "physSensitiveBridgeVolume" ) 
    return false;  
  */
  
  SensitiveDetectorHit* newHit = new SensitiveDetectorHit();

  newHit -> SetEdep(edep);
  
  // step length
  G4double len = aStep->GetStepLength();

  // only consider ions, protons, and neutrons for path length
  const G4ParticleDefinition* thisParticle = aStep ->
                      GetTrack() -> GetParticleDefinition() ;
  G4String particleName = thisParticle -> GetParticleName();
  G4int particleZ = thisParticle -> GetAtomicNumber();

  if ( (particleZ < 1) && (particleName != "neutron"))
    len = 0.;

  newHit->SetPath(len);
  
  fHitsCollection -> insert( newHit );

  return true;
}

void SensitiveDetector::EndOfEvent(G4HCofThisEvent*)
{
#ifdef ANALYSIS_USE
// Initialisation of total energy deposition per event to zero
 G4double totalEdepInOneEvent=0;
 G4double totalPathLengthInOneEvent=0;
 
 G4int NbHits = fHitsCollection->entries();
   //G4cout << "number of hits " <<NbHits << G4endl;
 
 G4double edep, len; 
 
 for (G4int i=0;i<NbHits;i++) 
 {
     edep = (*fHitsCollection)[i]->GetEdep();
     totalEdepInOneEvent = totalEdepInOneEvent + edep;
 
     len = (*fHitsCollection)[i]->GetPath();
     totalPathLengthInOneEvent = totalPathLengthInOneEvent + len;
 }
 
 G4int eventID = G4RunManager::GetRunManager() ->
     GetCurrentEvent() -> GetEventID();


 if (totalEdepInOneEvent!=0)
 {
     if (isMicrod==true) analysis -> StoreEnergyDeposition(totalEdepInOneEvent/keV, totalPathLengthInOneEvent/um, eventID);
     else analysis -> StoreSecondStageEnergyDeposition(totalEdepInOneEvent/keV, eventID); 
 }
#endif 
}
