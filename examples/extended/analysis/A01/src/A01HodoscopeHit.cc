// $Id: A01HodoscopeHit.cc,v 1.1 2002-11-13 07:23:40 duns Exp $
// --------------------------------------------------------------
// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
#include "A01HodoscopeHit.hh"
#include "G4VVisManager.hh"
#include "G4Circle.hh"
#include "G4Colour.hh"
#include "G4VisAttributes.hh"
#include "G4ios.hh"

G4Allocator<A01HodoscopeHit> A01HodoscopeHitAllocator;

A01HodoscopeHit::A01HodoscopeHit(G4int i,G4double t)
{
  id = i;
  time = t;
  pLogV = 0;
}

A01HodoscopeHit::~A01HodoscopeHit()
{;}

A01HodoscopeHit::A01HodoscopeHit(const A01HodoscopeHit &right)
{
  id = right.id;
  time = right.time;
  pos = right.pos;
  rot = right.rot;
  pLogV = right.pLogV;
}

const A01HodoscopeHit& A01HodoscopeHit::operator=(const A01HodoscopeHit &right)
{
  id = right.id;
  time = right.time;
  pos = right.pos;
  rot = right.rot;
  pLogV = right.pLogV;
  return *this;
}

int A01HodoscopeHit::operator==(const A01HodoscopeHit &right) const
{
  return 0;
}

void A01HodoscopeHit::Draw()
{
  G4VVisManager* pVVisManager = G4VVisManager::GetConcreteInstance();
  if(pVVisManager)
  {
    G4Transform3D trans(rot.inverse(),pos);
    G4VisAttributes attribs;
    const G4VisAttributes* pVA = pLogV->GetVisAttributes();
    if(pVA) attribs = *pVA;
    G4Colour colour(0.,1.,1.);
    attribs.SetColour(colour);
    attribs.SetForceWireframe(false);
    attribs.SetForceSolid(true);
    pVVisManager->Draw(*pLogV,attribs,trans);
  }
}

void A01HodoscopeHit::Print()
{
    G4cout << "  Hodoscope[" << id << "] " << time/ns << " (nsec)" << G4endl;
}



