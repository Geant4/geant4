// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4PDCofThisEvent.cc,v 1.3 1999/12/02 16:10:20 morita Exp $
// GEANT4 tag $Name: geant4-02-00 $
//

#include "G4PDCofThisEvent.hh"

G4PDCofThisEvent::G4PDCofThisEvent()
{
  DC.resize(10);
}

G4PDCofThisEvent::~G4PDCofThisEvent()
{;}

void G4PDCofThisEvent::AddDigitsCollection(G4int DCID,
                                           HepRef(G4PVDigitsCollection) aDC)
{
  if(DCID>=0)
  {
    if(DCID>=DC.size()) DC.resize(DCID+1);
    DC[DCID] = aDC;
  }
}


