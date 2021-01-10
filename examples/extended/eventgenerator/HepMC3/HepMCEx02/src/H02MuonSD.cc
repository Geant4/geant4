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
/// \file eventgenerator/HepMC/HepMCEx02/src/H02MuonSD.cc
/// \brief Implementation of the H02MuonSD class
//
//
#include "G4HCofThisEvent.hh"
#include "G4TouchableHistory.hh"
#include "G4Track.hh"
#include "G4Step.hh"
#include "G4VPhysicalVolume.hh"
#include "H02MuonHit.hh"
#include "H02MuonSD.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
H02MuonSD::H02MuonSD(G4String name)
  : G4VSensitiveDetector(name),
    fHitCollection(0)
{
  G4String HCname;
  collectionName.insert("muonHit");
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
H02MuonSD::~H02MuonSD()
{
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void H02MuonSD::Initialize(G4HCofThisEvent* HCE)
{
  static int HCID=-1;
  fHitCollection= new H02MuonHitsCollection(SensitiveDetectorName,
                                          collectionName[0]);
  if(HCID<0) HCID= GetCollectionID(0);
  HCE-> AddHitsCollection(HCID, fHitCollection);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
G4bool H02MuonSD::ProcessHits(G4Step* astep, G4TouchableHistory*)
{
  G4ParticleDefinition* particle= astep-> GetTrack()-> GetDefinition();
  if(particle-> GetPDGCharge() == 0.) return false;

  G4StepPoint* prestep= astep-> GetPreStepPoint();

  if(prestep-> GetStepStatus() != fGeomBoundary) return false;

  G4ThreeVector vmom= prestep-> GetMomentum();
  G4ThreeVector vpos= prestep-> GetPosition();
  G4double tof= prestep-> GetGlobalTime();

  G4VPhysicalVolume* volume= prestep-> GetPhysicalVolume();
  G4int id= volume-> GetCopyNo();
  if(volume-> GetName() == "ENDCAP_MUON_PV") id +=10;

  H02MuonHit* aHit=
    new H02MuonHit(id, particle-> GetParticleName(), vmom, vpos, tof);
  fHitCollection-> insert(aHit);
  return true;

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void H02MuonSD::EndOfEvent(G4HCofThisEvent*)
{
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void H02MuonSD::clear()
{
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void H02MuonSD::DrawAll()
{
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void H02MuonSD::PrintAll()
{
  G4int nHit= fHitCollection-> entries();
  G4cout << "------------------------------------------" << G4endl
         << "*** Muon System Hit (#hits=" << nHit << ")" << G4endl;
  fHitCollection-> PrintAllHits();
}
