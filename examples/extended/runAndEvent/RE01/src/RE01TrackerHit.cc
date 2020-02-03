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
/// \file runAndEvent/RE01/src/RE01TrackerHit.cc
/// \brief Implementation of the RE01TrackerHit class
//
//

#include "RE01TrackerHit.hh"

#include "G4VVisManager.hh"
#include "G4Circle.hh"
#include "G4Colour.hh"
#include "G4AttDefStore.hh"
#include "G4AttDef.hh"
#include "G4AttValue.hh"
#include "G4UIcommand.hh"
#include "G4UnitsTable.hh"
#include "G4VisAttributes.hh"
#include "G4SystemOfUnits.hh"    

G4ThreadLocal G4Allocator<RE01TrackerHit> * RE01TrackerHitAllocator = 0;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo...... 
RE01TrackerHit::RE01TrackerHit()
  : G4VHit(), fEdep(0.0), fPos(0),fTrackID(-1)
{;}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo...... 
RE01TrackerHit::~RE01TrackerHit()
{;}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo...... 
void RE01TrackerHit::Draw()
{
  G4VVisManager* pVVisManager = G4VVisManager::GetConcreteInstance();
  if(pVVisManager)
  {
    G4Circle circle(fPos);
    circle.SetScreenSize(0.04);
    circle.SetFillStyle(G4Circle::filled);
    G4Colour colour(1.,0.,0.);
    G4VisAttributes attribs(colour);
    circle.SetVisAttributes(attribs);
    pVVisManager->Draw(circle);
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo...... 
const std::map<G4String,G4AttDef>* RE01TrackerHit::GetAttDefs() const
{
  G4bool isNew;
  std::map<G4String,G4AttDef>* store
    = G4AttDefStore::GetInstance("RE01TrackerHit",isNew);
  if (isNew) {
    G4String hitType("HitType");
    (*store)[hitType] = G4AttDef(hitType,"Hit Type","Physics","","G4String");

    G4String trackID("TrackID");
    (*store)[trackID] = G4AttDef(trackID,"Track ID","Physics","","G4int");

    G4String energy("Energy");
    (*store)[energy] = G4AttDef(energy,"Energy Deposited","Physics",
                                "G4BestUnit","G4double");

    G4String eTrack("ETrack");
    (*store)[eTrack] = G4AttDef(eTrack,"Energy Deposited By Track","Physics",
                                "G4BestUnit","G4double");

    G4String pos("Pos");
    (*store)[pos] = G4AttDef(pos, "Position",
                      "Physics","G4BestUnit","G4ThreeVector");
  }
  return store;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo...... 
std::vector<G4AttValue>* RE01TrackerHit::CreateAttValues() const
{
  std::vector<G4AttValue>* values = new std::vector<G4AttValue>;

  values->push_back(G4AttValue("HitType","RE01TrackerHit",""));

  values->push_back
    (G4AttValue("TrackID",G4UIcommand::ConvertToString(fTrackID),""));

  values->push_back
    (G4AttValue("Energy",G4BestUnit(fEdep,"Energy"),""));

  G4double noEnergy = 0.*MeV;
  values->push_back
    (G4AttValue("ETrack",G4BestUnit(noEnergy,"Energy"),""));

  values->push_back
    (G4AttValue("Pos",G4BestUnit(fPos,"Length"),""));
  
   return values;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo...... 
void RE01TrackerHit::Print()
{
  G4cout << "TrackID " << fTrackID << "   Position " << fPos << "       : " 
         << fEdep/keV << " [keV]" << G4endl;
}
