// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: MyCalorimeterHit.cc,v 1.2 1999-08-27 10:19:07 johna Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//

#include "MyCalorimeterHit.hh"
#include "G4ios.hh"
#include "G4VVisManager.hh"
#include "G4Colour.hh"
#include "G4VisAttributes.hh"
#include "G4LogicalVolume.hh"

//G4Allocator<MyCalorimeterHit> MyCalorimeterHitAllocator;

MyCalorimeterHit::MyCalorimeterHit()
{pPhys=NULL;}

MyCalorimeterHit::MyCalorimeterHit(G4VPhysicalVolume* physVol)
:pPhys(physVol)
{;}

MyCalorimeterHit::~MyCalorimeterHit()
{;}

MyCalorimeterHit::MyCalorimeterHit(const MyCalorimeterHit &right)
{
  edep = right.edep;
  pos = right.pos;
  rot = right.rot;
  pPhys = right.pPhys;
}

const MyCalorimeterHit& MyCalorimeterHit::operator=(const MyCalorimeterHit &right)
{
  edep = right.edep;
  pos = right.pos;
  rot = right.rot;
  pPhys = right.pPhys;
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
    G4LogicalVolume* logVol = pPhys->GetLogicalVolume();
    const G4VisAttributes* pVA = logVol->GetVisAttributes();
    if(pVA) attribs = *pVA;
    G4Colour colour(1.,0.,0.);
    attribs.SetColour(colour);
    attribs.SetForceSolid(true);
    pVVisManager->Draw(*pPhys,attribs,trans);
  }
}

void MyCalorimeterHit::Print()
{
}


