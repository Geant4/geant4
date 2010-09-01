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

#include "Test2PhantomHit.hh"
#include "G4LogicalVolume.hh"
#include "G4ThreeVector.hh"
#include "G4Transform3D.hh"
#include "G4RotationMatrix.hh"
#include "G4UIcommand.hh"
#include "G4UnitsTable.hh"

#include "G4VVisManager.hh"
#include "G4VisAttributes.hh"
#include "G4AttValue.hh"
#include "G4AttDef.hh"
#include "G4AttCheck.hh"
#include "G4Colour.hh"

G4Allocator<Test2PhantomHit> Test2PhantomHitAllocator;

Test2PhantomHit::Test2PhantomHit() {
  fpLogV = 0;
}

Test2PhantomHit::Test2PhantomHit(G4LogicalVolume * logVol,
				 G4int & x, G4int & y, G4int & z)
  : fXCellID(x), fYCellID(y), fZCellID(z), fpLogV(logVol) {
  ;
}

Test2PhantomHit::~Test2PhantomHit() {
  ;
}

Test2PhantomHit::Test2PhantomHit(const Test2PhantomHit & right)
  : G4VHit() {

  fXCellID = right.fXCellID;
  fYCellID = right.fYCellID;
  fZCellID = right.fZCellID;
  fEdep = right.fEdep;
  fTrackLength = right.fTrackLength;
  fParticleName = right.fParticleName;
  fPos = right.fPos;
  fRotMat = right.fRotMat;
  fpLogV = right.fpLogV;
}

const Test2PhantomHit & Test2PhantomHit::operator=(const Test2PhantomHit & right) {

  fXCellID = right.fXCellID;
  fYCellID = right.fYCellID;
  fZCellID = right.fZCellID;
  fEdep = right.fEdep;
  fTrackLength = right.fTrackLength;
  fParticleName = right.fParticleName;
  fPos = right.fPos;
  fRotMat = right.fRotMat;
  fpLogV = right.fpLogV;

  return *this;
}

G4int Test2PhantomHit::operator==(const Test2PhantomHit &right) const {

  return ((fXCellID == right.fXCellID) &&
	  (fYCellID == right.fYCellID) &&
	  (fZCellID == right.fZCellID));
}

std::map<G4String, G4AttDef> Test2PhantomHit::fAttDefs;

void Test2PhantomHit::Draw() {

  G4VVisManager* pVVisManager = G4VVisManager::GetConcreteInstance();
  if(pVVisManager) {
    G4Transform3D trans(fRotMat, fPos);
    G4VisAttributes attribs;
    const G4VisAttributes* pVA = fpLogV->GetVisAttributes();
    if(pVA) attribs = *pVA;
    G4Colour colour(1.,0.,0.);
    attribs.SetColour(colour);
    attribs.SetForceSolid(true);
    pVVisManager->Draw(*fpLogV, attribs, trans);
  }
}

const std::map<G4String,G4AttDef>* Test2PhantomHit::GetAttDefs() const {

  // G4AttDefs have to have long life.  Use static member...
  if (fAttDefs.empty()) {
    fAttDefs["HitType"] =
      G4AttDef("HitType","Type of hit","Physics","","G4String");
    fAttDefs["XID"] = G4AttDef("XID","X Cell ID","Physics","","G4int");
    fAttDefs["YID"] = G4AttDef("YID","Y Cell ID","Physics","","G4int");
    fAttDefs["ZID"] = G4AttDef("ZID","Z Cell ID","Physics","","G4int");
    fAttDefs["EDep"] =
      G4AttDef("EDep","Energy deposited","Physics","G4BestUnit","G4double");
  }
  return &fAttDefs;
}

std::vector<G4AttValue>* Test2PhantomHit::CreateAttValues() const {

  // Create expendable G4AttsValues for picking...
  std::vector<G4AttValue>* attValues = new std::vector<G4AttValue>;
  attValues->push_back(G4AttValue("HitType","Test2PhantomHit",""));
  attValues->push_back(G4AttValue("XID",G4UIcommand::ConvertToString(fXCellID),""));
  attValues->push_back(G4AttValue("YID",G4UIcommand::ConvertToString(fYCellID),""));
  attValues->push_back(G4AttValue("ZID",G4UIcommand::ConvertToString(fZCellID),""));
  attValues->push_back(G4AttValue("EDep",G4BestUnit(fEdep,"Energy"),""));
  //G4cout << "Checking...\n" << G4AttCheck(attValues, GetAttDefs());
  return attValues;
}

void Test2PhantomHit::Print()
{;}


