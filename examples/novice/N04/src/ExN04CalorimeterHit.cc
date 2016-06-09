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

#include "ExN04CalorimeterHit.hh"
#include "G4ios.hh"
#include "G4VVisManager.hh"
#include "G4Colour.hh"
#include "G4VisAttributes.hh"
#include "G4LogicalVolume.hh"
#include "G4UIcommand.hh"
#include "G4UnitsTable.hh"
#include "G4AttValue.hh"
#include "G4AttDef.hh"
#include "G4AttCheck.hh"

G4Allocator<ExN04CalorimeterHit> ExN04CalorimeterHitAllocator;

ExN04CalorimeterHit::ExN04CalorimeterHit()
{pLogV=0;}

ExN04CalorimeterHit::ExN04CalorimeterHit(G4LogicalVolume* logVol,G4int z,G4int phi)
: ZCellID(z), PhiCellID(phi), pLogV(logVol)
{;}

ExN04CalorimeterHit::~ExN04CalorimeterHit()
{;}

ExN04CalorimeterHit::ExN04CalorimeterHit(const ExN04CalorimeterHit &right)
  : G4VHit()
{
  ZCellID = right.ZCellID;
  PhiCellID = right.PhiCellID;
  edep = right.edep;
  pos = right.pos;
  rot = right.rot;
  pLogV = right.pLogV;
}

const ExN04CalorimeterHit& ExN04CalorimeterHit::operator=(const ExN04CalorimeterHit &right)
{
  ZCellID = right.ZCellID;
  PhiCellID = right.PhiCellID;
  edep = right.edep;
  pos = right.pos;
  rot = right.rot;
  pLogV = right.pLogV;
  return *this;
}

G4int ExN04CalorimeterHit::operator==(const ExN04CalorimeterHit &right) const
{
  return ((ZCellID==right.ZCellID)&&(PhiCellID==right.PhiCellID));
}

std::map<G4String,G4AttDef> ExN04CalorimeterHit::fAttDefs;

void ExN04CalorimeterHit::Draw()
{
  G4VVisManager* pVVisManager = G4VVisManager::GetConcreteInstance();
  if(pVVisManager)
  {
    G4Transform3D trans(rot,pos);
    G4VisAttributes attribs;
    const G4VisAttributes* pVA = pLogV->GetVisAttributes();
    if(pVA) attribs = *pVA;
    G4Colour colour(1.,0.,0.);
    attribs.SetColour(colour);
    attribs.SetForceSolid(true);
    pVVisManager->Draw(*pLogV,attribs,trans);
  }
}

const std::map<G4String,G4AttDef>* ExN04CalorimeterHit::GetAttDefs() const
{
  // G4AttDefs have to have long life.  Use static member...
  if (fAttDefs.empty()) {
    fAttDefs["HitType"] =
      G4AttDef("HitType","Type of hit","Physics","","G4String");
    fAttDefs["ZID"] = G4AttDef("ZID","Z Cell ID","Physics","","G4int");
    fAttDefs["PhiID"] = G4AttDef("PhiID","Phi Cell ID","Physics","","G4int");
    fAttDefs["EDep"] =
      G4AttDef("EDep","Energy deposited","Physics","G4BestUnit","G4double");
  }
  return &fAttDefs;
}

std::vector<G4AttValue>* ExN04CalorimeterHit::CreateAttValues() const
{
  // Create expendable G4AttsValues for picking...
  std::vector<G4AttValue>* attValues = new std::vector<G4AttValue>;
  attValues->push_back
    (G4AttValue("HitType","ExN04CalorimeterHit",""));
  attValues->push_back
    (G4AttValue("ZID",G4UIcommand::ConvertToString(ZCellID),""));
  attValues->push_back
    (G4AttValue("PhiID",G4UIcommand::ConvertToString(PhiCellID),""));
  attValues->push_back
    (G4AttValue("EDep",G4BestUnit(edep,"Energy"),""));
  //G4cout << "Checking...\n" << G4AttCheck(attValues, GetAttDefs());
  return attValues;
}

void ExN04CalorimeterHit::Print()
{;}


