// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4THitsCollection.cc,v 1.1.10.1 1999/12/07 20:47:48 gunter Exp $
// GEANT4 tag $Name: geant4-02-00 $
//

#include "G4THitsCollection.hh"

G4Allocator<G4HitsCollection> anHCAllocator;

G4HitsCollection::G4HitsCollection()
{;}

G4HitsCollection::G4HitsCollection(G4String detName,G4String colNam)
: G4VHitsCollection(detName,colNam)
{;}

G4HitsCollection::~G4HitsCollection()
{;}

int G4HitsCollection::operator==(const G4HitsCollection &right) const
{ return (collectionName==right.collectionName); }

