// $Id: A01EmCalorimeterHit.cc,v 1.1 2002-11-13 07:23:08 duns Exp $
// --------------------------------------------------------------
// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//

#include "A01EmCalorimeterHit.hh"
#include "G4VVisManager.hh"
#include "G4Colour.hh"
#include "G4VisAttributes.hh"
#include "G4LogicalVolume.hh"
#include "G4ios.hh"

G4Allocator<A01EmCalorimeterHit> A01EmCalorimeterHitAllocator;

A01EmCalorimeterHit::A01EmCalorimeterHit()
{
  cellID = -1;
  edep = 0.;
  pLogV = 0;
}

A01EmCalorimeterHit::A01EmCalorimeterHit(G4int z)
{
  cellID = z;
  edep = 0.;
  pLogV = 0;
}

A01EmCalorimeterHit::~A01EmCalorimeterHit()
{;}

A01EmCalorimeterHit::A01EmCalorimeterHit(const A01EmCalorimeterHit &right)
{
  cellID = right.cellID;
  edep = right.edep;
  pos = right.pos;
  rot = right.rot;
  pLogV = right.pLogV;
}

const A01EmCalorimeterHit& A01EmCalorimeterHit::operator=(const A01EmCalorimeterHit &right)
{
  cellID = right.cellID;
  edep = right.edep;
  pos = right.pos;
  rot = right.rot;
  pLogV = right.pLogV;
  return *this;
}

int A01EmCalorimeterHit::operator==(const A01EmCalorimeterHit &right) const
{
  return (cellID==right.cellID);
}

void A01EmCalorimeterHit::Draw()
{
  G4VVisManager* pVVisManager = G4VVisManager::GetConcreteInstance();
  if(pVVisManager&&(edep>0.))
  {
    // Draw a calorimeter cell with a color corresponding to its energy deposit
    G4Transform3D trans(rot.inverse(),pos);
    G4VisAttributes attribs;
    const G4VisAttributes* pVA = pLogV->GetVisAttributes();
    if(pVA) attribs = *pVA;
    G4double rcol = edep/(0.7*GeV);
    if(rcol>1.) rcol = 1.;
    G4Colour colour(rcol,0.,0.);
    attribs.SetColour(colour);
    attribs.SetForceWireframe(false);
    attribs.SetForceSolid(true);
    pVVisManager->Draw(*pLogV,attribs,trans);
  }
}

void A01EmCalorimeterHit::Print()
{
  G4cout << "  Cell[" << cellID << "] " << edep/MeV << " (MeV)" << G4endl;
}


