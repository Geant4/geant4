// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4DCofThisEvent.cc,v 1.4 2001/02/08 06:07:14 asaim Exp $
// GEANT4 tag $Name: geant4-03-01 $
//

#include "G4DCofThisEvent.hh"

G4Allocator<G4DCofThisEvent> anDCoTHAllocator;

G4DCofThisEvent::G4DCofThisEvent()
{
  DC = new G4std::vector<G4VDigiCollection*>;
}

G4DCofThisEvent::G4DCofThisEvent(G4int cap)
{
  DC = new G4std::vector<G4VDigiCollection*>;
  for(int i=0;i<cap;i++)
  {
    DC->push_back((G4VDigiCollection*)NULL);
  }
}

G4DCofThisEvent::~G4DCofThisEvent()
{
  //DC->clearAndDestroy();
  for(G4int i=0;i<DC->size();i++)
  { delete (*DC)[i]; }
  DC->clear();
  delete DC;
}

void G4DCofThisEvent::AddDigiCollection(G4int DCID,G4VDigiCollection * aDC)
{
  if(DCID>=0 && DCID<DC->size())
  { (*DC)[DCID] = aDC; }
}

