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
// ====================================================================
//
//   H02MuonSD.cc
//   $Id: H02MuonSD.cc,v 1.5 2006-07-05 12:04:13 gcosmo Exp $
//
// ====================================================================

#include "H02MuonSD.hh"
#include "H02MuonHit.hh"
#include "G4HCofThisEvent.hh"
#include "G4TouchableHistory.hh"
#include "G4VPhysicalVolume.hh"
#include "G4Track.hh"
#include "G4Step.hh"
#include "G4ios.hh"

/////////////////////////////////
H02MuonSD::H02MuonSD(G4String name) 
  : G4VSensitiveDetector(name)
/////////////////////////////////
{
  G4String HCname;
  collectionName.insert("muonHit");
}

///////////////////////////////////////////////
H02MuonSD::~H02MuonSD()
///////////////////////////////////////////////
{
}

///////////////////////////////////////////////
void H02MuonSD::Initialize(G4HCofThisEvent* HCE)
///////////////////////////////////////////////
{
  static int HCID=-1;
  hitCollection= new H02MuonHitsCollection(SensitiveDetectorName, 
					  collectionName[0]); 
  if(HCID<0) HCID= GetCollectionID(0);
  HCE-> AddHitsCollection(HCID, hitCollection);
}

///////////////////////////////////////////////////////////////////////
G4bool H02MuonSD::ProcessHits(G4Step* astep, G4TouchableHistory*)
///////////////////////////////////////////////////////////////////////
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
  hitCollection-> insert(aHit);
  return true;

}

///////////////////////////////////////////////
void H02MuonSD::EndOfEvent(G4HCofThisEvent*)
///////////////////////////////////////////////
{
}

//////////////////////
void H02MuonSD::clear()
//////////////////////
{
} 

////////////////////////
void H02MuonSD::DrawAll()
////////////////////////
{
} 

/////////////////////////
void H02MuonSD::PrintAll()
/////////////////////////
{
  G4int nHit= hitCollection-> entries();
  G4cout << "------------------------------------------" << G4endl
         << "*** Muon System Hit (#hits=" << nHit << ")" << G4endl;
  hitCollection-> PrintAllHits();
} 

