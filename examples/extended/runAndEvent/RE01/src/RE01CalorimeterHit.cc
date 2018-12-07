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
/// \file runAndEvent/RE01/src/RE01CalorimeterHit.cc
/// \brief Implementation of the RE01CalorimeterHit class
//
//
//

#include "RE01CalorimeterHit.hh"
#include "G4ios.hh"
#include "G4VVisManager.hh"
#include "G4Colour.hh"
#include "G4AttDefStore.hh"
#include "G4AttDef.hh"
#include "G4AttValue.hh"
#include "G4UIcommand.hh"
#include "G4UnitsTable.hh"
#include "G4VisAttributes.hh"
#include "G4LogicalVolume.hh"
#include "G4SystemOfUnits.hh"    

G4ThreadLocal G4Allocator<RE01CalorimeterHit> * RE01CalorimeterHitAllocator = 0;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
RE01CalorimeterHit::RE01CalorimeterHit(G4LogicalVolume* logVol,
                                       G4int z, G4int phi)
  :G4VHit(), fZCellID(z), fPhiCellID(phi), fEdep(0.0), 
   fPos(0),fRot(0.,0.,0.),fPLogV(logVol),fEdepByATrack(0.0),fTrackInfo()
{;}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
RE01CalorimeterHit::~RE01CalorimeterHit()
{;}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void RE01CalorimeterHit::Draw()
{
  G4VVisManager* pVVisManager = G4VVisManager::GetConcreteInstance();
  if(pVVisManager)
  {
    G4Transform3D trans(fRot,fPos);
    G4VisAttributes attribs;
    const G4VisAttributes* pVA = fPLogV->GetVisAttributes();
    if(pVA) attribs = *pVA;
    G4Colour colour(1.,0.,0.);
    attribs.SetColour(colour);
    attribs.SetForceWireframe(false);
    attribs.SetForceSolid(true);
    pVVisManager->Draw(*fPLogV,attribs,trans);
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
const std::map<G4String,G4AttDef>* RE01CalorimeterHit::GetAttDefs() const
{
  G4bool isNew;
  std::map<G4String,G4AttDef>* store
    = G4AttDefStore::GetInstance("RE01CalorimeterHit",isNew);
  if (isNew) {
    G4String hitType("HitType");
    (*store)[hitType] = G4AttDef(hitType,"Hit Type","Physics","","G4String");

    G4String zCellID("ZCellID");
    (*store)[zCellID] = G4AttDef(zCellID,"Z Cell ID","Physics","","G4int");

    G4String phiCellID("PhiCellID");
    (*store)[phiCellID] = G4AttDef(phiCellID,"Phi Cell ID","Physics","","G4int");

    G4String energy("Energy");
    (*store)[energy] = G4AttDef(energy,"Energy Deposited","Physics","G4BestUnit",
                                "G4double");

    G4String eTrack("ETrack");
    (*store)[eTrack] = G4AttDef(eTrack,"Energy Deposited By Track","Physics",
                                "G4BestUnit","G4double");

    G4String pos("Pos");
    (*store)[pos] = G4AttDef(pos, "Position",
                      "Physics","G4BestUnit","G4ThreeVector");

    G4String lvol("LVol");
    (*store)[lvol] = G4AttDef(lvol,"Logical Volume","Physics","","G4String");
  }

  return store;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
std::vector<G4AttValue>* RE01CalorimeterHit::CreateAttValues() const
{
  std::vector<G4AttValue>* values = new std::vector<G4AttValue>;

  values->push_back(G4AttValue("HitType","RE01CalorimeterHit",""));

  values->push_back
    (G4AttValue("ZCellID",G4UIcommand::ConvertToString(fZCellID),""));

  values->push_back
    (G4AttValue("PhiCellID",G4UIcommand::ConvertToString(fPhiCellID),""));

  values->push_back
    (G4AttValue("Energy",G4BestUnit(fEdep,"Energy"),""));

  values->push_back
    (G4AttValue("ETrack",G4BestUnit(fEdepByATrack,"Energy"),""));

  values->push_back
    (G4AttValue("Pos",G4BestUnit(fPos,"Length"),""));

  if (fPLogV)
    values->push_back
      (G4AttValue("LVol",fPLogV->GetName(),""));
  else
    values->push_back
      (G4AttValue("LVol"," ",""));

  return values;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void RE01CalorimeterHit::Print()
{
  G4cout << "Cell[" << fZCellID << "," << fPhiCellID << "]    " 
         << fEdep/GeV << " [GeV]" << G4endl;
}


