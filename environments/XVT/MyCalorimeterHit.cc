// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: MyCalorimeterHit.cc,v 1.1 1999-01-07 16:04:58 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//

#include "MyCalorimeterHit.hh"
#include "G4ios.hh"
#include "G4VVisManager.hh"
#include "G4Colour.hh"
#include "G4VisAttributes.hh"
#include "G4LogicalVolume.hh"

G4Allocator<MyCalorimeterHit> MyCalorimeterHitAllocator;

MyCalorimeterHit::MyCalorimeterHit()
{pLogV=NULL;}

MyCalorimeterHit::MyCalorimeterHit(G4LogicalVolume* logVol)
:pLogV(logVol)
{;}

MyCalorimeterHit::~MyCalorimeterHit()
{;}

MyCalorimeterHit::MyCalorimeterHit(const MyCalorimeterHit &right)
{
  edep = right.edep;
  pos = right.pos;
  rot = right.rot;
  pLogV = right.pLogV;
}

const MyCalorimeterHit& MyCalorimeterHit::operator=(const MyCalorimeterHit &right)
{
  edep = right.edep;
  pos = right.pos;
  rot = right.rot;
  pLogV = right.pLogV;
  return *this;
}

int MyCalorimeterHit::operator==(const MyCalorimeterHit &right) const
{
  return 0;
}

void MyCalorimeterHit::Draw()
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

void MyCalorimeterHit::Print()
{
}


