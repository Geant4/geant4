//
// ********************************************************************
// * DISCLAIMER                                                       *
// *                                                                  *
// * The following disclaimer summarizes all the specific disclaimers *
// * of contributors to this software. The specific disclaimers,which *
// * govern, are listed with their locations in:                      *
// *   http://cern.ch/geant4/license                                  *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.                                                             *
// *                                                                  *
// * This  code  implementation is the  intellectual property  of the *
// * GEANT4 collaboration.                                            *
// * By copying,  distributing  or modifying the Program (or any work *
// * based  on  the Program)  you indicate  your  acceptance of  this *
// * statement, and all its terms.                                    *
// ********************************************************************

// ====================================================================
//
//   XXMuonSD.cc
//   $Id: XXMuonSD.cc,v 1.1 2002-04-29 20:45:01 asaim Exp $
//
// ====================================================================

#include "XXMuonSD.hh"
#include "XXMuonHit.hh"
#include "G4HCofThisEvent.hh"
#include "G4TouchableHistory.hh"
#include "G4VPhysicalVolume.hh"
#include "G4Track.hh"
#include "G4Step.hh"
#include "G4ios.hh"

/////////////////////////////////
XXMuonSD::XXMuonSD(G4String name) 
  : G4VSensitiveDetector(name)
/////////////////////////////////
{
  G4String HCname;
  collectionName.insert("muonHit");
}

///////////////////////////////////////////////
XXMuonSD::~XXMuonSD()
///////////////////////////////////////////////
{
}

///////////////////////////////////////////////
void XXMuonSD::Initialize(G4HCofThisEvent* HCE)
///////////////////////////////////////////////
{
  static int HCID=-1;
  hitCollection= new XXMuonHitsCollection(SensitiveDetectorName, 
					  collectionName[0]); 
  if(HCID<0) HCID= GetCollectionID(0);
  HCE-> AddHitsCollection(HCID, hitCollection);
}

///////////////////////////////////////////////////////////////////////
G4bool XXMuonSD::ProcessHits(G4Step* astep, G4TouchableHistory* ROhist)
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
  
  XXMuonHit* aHit= 
    new XXMuonHit(id, particle-> GetParticleName(), vmom, vpos, tof);  
  hitCollection-> insert(aHit);
  return true;

}

///////////////////////////////////////////////
void XXMuonSD::EndOfEvent(G4HCofThisEvent* HCE)
///////////////////////////////////////////////
{
}

//////////////////////
void XXMuonSD::clear()
//////////////////////
{
} 

////////////////////////
void XXMuonSD::DrawAll()
////////////////////////
{
  hitCollection-> DrawAllHits();
} 

/////////////////////////
void XXMuonSD::PrintAll()
/////////////////////////
{
  G4int nHit= hitCollection-> entries();
  G4cout << "------------------------------------------" << G4endl
         << "*** Muon System Hit (#hits=" << nHit << ")" << G4endl;
  hitCollection-> PrintAllHits();
} 

