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
//
//
// $Id: RE01CalorimeterHit.cc,v 1.3 2005/06/07 10:53:17 perl Exp $
// GEANT4 tag $Name: geant4-08-00 $
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

G4Allocator<RE01CalorimeterHit> RE01CalorimeterHitAllocator;

RE01CalorimeterHit::RE01CalorimeterHit()
{pLogV=0;}

RE01CalorimeterHit::RE01CalorimeterHit(G4LogicalVolume* logVol,G4int z,G4int phi)
: ZCellID(z), PhiCellID(phi), pLogV(logVol)
{;}

RE01CalorimeterHit::~RE01CalorimeterHit()
{;}

RE01CalorimeterHit::RE01CalorimeterHit(const RE01CalorimeterHit &right)
  : G4VHit()
{
  ZCellID = right.ZCellID;
  PhiCellID = right.PhiCellID;
  edep = right.edep;
  pos = right.pos;
  rot = right.rot;
  pLogV = right.pLogV;
  edepByATrack = right.edepByATrack;
  trackInfo = right.trackInfo;
}

const RE01CalorimeterHit& RE01CalorimeterHit::operator=(const RE01CalorimeterHit &right)
{
  ZCellID = right.ZCellID;
  PhiCellID = right.PhiCellID;
  edep = right.edep;
  pos = right.pos;
  rot = right.rot;
  pLogV = right.pLogV;
  edepByATrack = right.edepByATrack;
  trackInfo = right.trackInfo;
  return *this;
}

G4int RE01CalorimeterHit::operator==(const RE01CalorimeterHit &right) const
{
  return ((ZCellID==right.ZCellID)&&(PhiCellID==right.PhiCellID));
}

void RE01CalorimeterHit::Draw()
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
    attribs.SetForceWireframe(false);
    attribs.SetForceSolid(true);
    pVVisManager->Draw(*pLogV,attribs,trans);
  }
}

const std::map<G4String,G4AttDef>* RE01CalorimeterHit::GetAttDefs() const
{
  G4bool isNew;
  std::map<G4String,G4AttDef>* store
    = G4AttDefStore::GetInstance("RE01CalorimeterHit",isNew);
  if (isNew) {
    G4String HitType("HitType");
    (*store)[HitType] = G4AttDef(HitType,"Hit Type","Bookkeeping","","G4String");

    G4String TrackID("TrackID");
    (*store)[TrackID] = G4AttDef(TrackID,"Track ID","Bookkeeping","","G4int");

    G4String ZCellID("ZCellID");
    (*store)[ZCellID] = G4AttDef(ZCellID,"Z Cell ID","Bookkeeping","","G4int");

    G4String PhiCellID("PhiCellID");
    (*store)[PhiCellID] = G4AttDef(PhiCellID,"Phi Cell ID","Bookkeeping","","G4int");

    G4String Energy("Energy");
    (*store)[Energy] = G4AttDef(Energy,"Energy Deposited","Physics","G4BestUnit","G4double");

    G4String ETrack("ETrack");
    (*store)[ETrack] = G4AttDef(ETrack,"Energy Deposited By Track","Physics","G4BestUnit","G4double");

    G4String Pos("Pos");
    (*store)[Pos] = G4AttDef(Pos, "Position",
		      "Physics","G4BestUnit","G4ThreeVector");

    G4String LVol("LVol");
    (*store)[LVol] = G4AttDef(LVol,"Logical Volume","Bookkeeping","","G4String");
  }
  return store;
}

std::vector<G4AttValue>* RE01CalorimeterHit::CreateAttValues() const
{
  std::vector<G4AttValue>* values = new std::vector<G4AttValue>;

  values->push_back(G4AttValue("HitType","RE01CalorimeterHit",""));

  values->push_back
    (G4AttValue("TrackID"," ",""));

  values->push_back
    (G4AttValue("ZCellID",G4UIcommand::ConvertToString(ZCellID),""));

  values->push_back
    (G4AttValue("PhiCellID",G4UIcommand::ConvertToString(PhiCellID),""));

  values->push_back
    (G4AttValue("Energy",G4BestUnit(edep,"Energy"),""));

  values->push_back
    (G4AttValue("ETrack",G4BestUnit(edepByATrack,"Energy"),""));

  values->push_back
    (G4AttValue("Pos",G4BestUnit(pos,"Length"),""));

  if (pLogV)
    values->push_back
      (G4AttValue("LVol",pLogV->GetName(),""));
  else
    values->push_back
      (G4AttValue("LVol"," ",""));

  return values;
}

void RE01CalorimeterHit::Print()
{
  G4cout << "Cell[" << ZCellID << "," << PhiCellID << "]    " << edep/GeV << " [GeV]" << G4endl;
}


