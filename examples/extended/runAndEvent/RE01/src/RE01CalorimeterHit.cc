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
// $Id: RE01CalorimeterHit.cc,v 1.1 2004/11/26 07:37:41 asaim Exp $
// GEANT4 tag $Name: geant4-07-00-cand-01 $
//


#include "RE01CalorimeterHit.hh"
#include "G4ios.hh"
#include "G4VVisManager.hh"
#include "G4Colour.hh"
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

void RE01CalorimeterHit::Print()
{
  G4cout << "Cell[" << ZCellID << "," << PhiCellID << "]    " << edep/GeV << " [GeV]" << G4endl;
}


