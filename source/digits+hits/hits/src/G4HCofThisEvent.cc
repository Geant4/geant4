// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4HCofThisEvent.cc,v 1.1 1999/01/07 16:06:31 gunter Exp $
// GEANT4 tag $Name: geant4-00-01 $
//

#include "G4HCofThisEvent.hh"

G4Allocator<G4HCofThisEvent> anHCoTHAllocator;

G4HCofThisEvent::G4HCofThisEvent()
{
  HC = new RWTPtrOrderedVector<G4VHitsCollection>;
}

G4HCofThisEvent::G4HCofThisEvent(G4int cap)
{
  HC = new RWTPtrOrderedVector<G4VHitsCollection>;
  for(int i=0;i<cap;i++)
  {
    HC->insert((G4VHitsCollection*)NULL);
  }
}

G4HCofThisEvent::~G4HCofThisEvent()
{
  HC->clearAndDestroy();
  delete HC;
}

void G4HCofThisEvent::AddHitsCollection(G4int HCID,G4VHitsCollection * aHC)
{
  if(HCID>=0 && HCID<HC->entries())
  { (*HC)[HCID] = aHC; }
}


