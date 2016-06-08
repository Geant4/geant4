// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4VHitsCollection.cc,v 1.1.10.1 1999/12/07 20:47:49 gunter Exp $
// GEANT4 tag $Name: geant4-01-01 $
//

// G4VHitsCollection

#include "G4VHitsCollection.hh"

G4VHitsCollection::G4VHitsCollection()
{
  collectionName = "Unknown";
  SDname = "Unknown";
}

G4VHitsCollection::G4VHitsCollection(G4String detName,G4String colNam)
{
  collectionName = colNam;
  SDname = detName;
}

G4VHitsCollection::~G4VHitsCollection()
{ ; }

int G4VHitsCollection::operator==(const G4VHitsCollection &right) const
{ 
  return ((collectionName==right.collectionName)
        &&(SDname==right.SDname));
}

void G4VHitsCollection::DrawAllHits()
{;}

void G4VHitsCollection::PrintAllHits()
{;}

