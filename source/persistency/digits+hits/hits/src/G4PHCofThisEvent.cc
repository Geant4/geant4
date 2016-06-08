// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4PHCofThisEvent.cc,v 1.6 1999/12/02 16:10:22 morita Exp $
// GEANT4 tag $Name: geant4-01-01 $
//

#include "G4PHCofThisEvent.hh"

G4PHCofThisEvent::G4PHCofThisEvent()
{
  HC.resize(10);
}

G4PHCofThisEvent::~G4PHCofThisEvent()
{;}

void G4PHCofThisEvent::AddHitsCollection(G4int HCID,
                                         HepRef(G4PVHitsCollection) aHC)
{
  if(HCID>=0)
  {
    if(HCID>=HC.size()) HC.resize(HCID+1);
    HC[HCID] = aHC;
  }
}

