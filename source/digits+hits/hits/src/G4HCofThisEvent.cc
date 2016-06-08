// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4HCofThisEvent.cc,v 1.4 2001/02/08 06:07:15 asaim Exp $
// GEANT4 tag $Name: geant4-03-01 $
//

#include "G4HCofThisEvent.hh"

G4Allocator<G4HCofThisEvent> anHCoTHAllocator;

G4HCofThisEvent::G4HCofThisEvent()
{
  HC = new G4std::vector<G4VHitsCollection*>;
}

G4HCofThisEvent::G4HCofThisEvent(G4int cap)
{
  HC = new G4std::vector<G4VHitsCollection*>;
  for(int i=0;i<cap;i++)
  {
    HC->push_back((G4VHitsCollection*)NULL);
  }
}

G4HCofThisEvent::~G4HCofThisEvent()
{
  //HC->clearAndDestroy();
  for(G4int i=0;i<HC->size();i++)
  { delete (*HC)[i]; }
  HC->clear();
  delete HC;
}

void G4HCofThisEvent::AddHitsCollection(G4int HCID,G4VHitsCollection * aHC)
{
  if(HCID>=0 && HCID<HC->size())
  { (*HC)[HCID] = aHC; }
}


