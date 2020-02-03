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
/// \file eventgenerator/HepMC/HepMCEx01/src/ExN04MuonSD.cc
/// \brief Implementation of the ExN04MuonSD class
//
//

#include "G4HCofThisEvent.hh"
#include "G4ios.hh"
#include "G4Step.hh"
#include "G4SystemOfUnits.hh"
#include "G4TouchableHistory.hh"
#include "G4Track.hh"
#include "ExN04MuonSD.hh"
#include "ExN04MuonHit.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
ExN04MuonSD::ExN04MuonSD(G4String name)
  : G4VSensitiveDetector(name), 
    fMuonCollection(NULL), fPositionResolution(5*cm)
{
  G4String HCname;
  collectionName.insert(HCname="muonCollection");
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
ExN04MuonSD::~ExN04MuonSD(){;}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void ExN04MuonSD::Initialize(G4HCofThisEvent*HCE)
{
  static int HCID = -1;
  fMuonCollection = new ExN04MuonHitsCollection
                   (SensitiveDetectorName,collectionName[0]);
  if(HCID<0)
  { HCID = GetCollectionID(0); }
  HCE->AddHitsCollection(HCID,fMuonCollection);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
G4bool ExN04MuonSD::ProcessHits(G4Step*aStep,G4TouchableHistory*)
{
  G4double edep = aStep->GetTotalEnergyDeposit();

  if(edep==0.) return true;

  ExN04MuonHit* aHit;
  int nHit = fMuonCollection->entries();
  G4ThreeVector hitpos = aStep->GetPreStepPoint()->GetPosition();
  for(int i=0;i<nHit;i++) {
    aHit = (*fMuonCollection)[i];
    G4ThreeVector pos = aHit->GetPos();
    G4double dist2 = sqr(pos.x()-hitpos.x())
                    +sqr(pos.y()-hitpos.y())+sqr(pos.z()-hitpos.z());
    if(dist2<=sqr(fPositionResolution))
    aHit->AddEdep(edep);
    return true;
  }

  aHit = new ExN04MuonHit();
  aHit->SetEdep( edep );
  aHit->SetPos( aStep->GetPreStepPoint()->GetPosition() );
  fMuonCollection->insert( aHit );

  return true;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void ExN04MuonSD::EndOfEvent(G4HCofThisEvent*)
{
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void ExN04MuonSD::clear()
{
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void ExN04MuonSD::DrawAll()
{
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void ExN04MuonSD::PrintAll()
{
}
