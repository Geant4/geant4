// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4DCofThisEvent.cc,v 1.1 1999/01/07 16:06:29 gunter Exp $
// GEANT4 tag $Name: geant4-00-01 $
//

#include "G4DCofThisEvent.hh"

G4Allocator<G4DCofThisEvent> anDCoTHAllocator;

G4DCofThisEvent::G4DCofThisEvent()
{
  DC = new RWTPtrOrderedVector<G4VDigiCollection>;
}

G4DCofThisEvent::G4DCofThisEvent(G4int cap)
{
  DC = new RWTPtrOrderedVector<G4VDigiCollection>;
  for(int i=0;i<cap;i++)
  {
    DC->insert((G4VDigiCollection*)NULL);
  }
}

G4DCofThisEvent::~G4DCofThisEvent()
{
  DC->clearAndDestroy();
  delete DC;
}

void G4DCofThisEvent::AddDigiCollection(G4int DCID,G4VDigiCollection * aDC)
{
  if(DCID>=0 && DCID<DC->entries())
  { (*DC)[DCID] = aDC; }
}

