// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: MyCalorimeterHitsCollection.cc,v 1.2 1999-12-15 14:48:45 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//

#include "MyCalorimeterHitsCollection.hh"

MyCalorimeterHitsCollection::MyCalorimeterHitsCollection()
{}

MyCalorimeterHitsCollection::MyCalorimeterHitsCollection(G4String aName,
                               G4VSensitiveDetector * theSD)
:G4VHitsCollection(aName,theSD)
{}

MyCalorimeterHitsCollection::~MyCalorimeterHitsCollection()
{}

void MyCalorimeterHitsCollection::DrawAllHits()
{
  G4int n_hit = theCollection.entries();
  for(G4int i=0;i<n_hit;i++)
  { theCollection[i].Draw(); }
}

void MyCalorimeterHitsCollection::PrintAllHits()
{
  G4int n_hit = theCollection.entries();
  for(G4int i=0;i<n_hit;i++)
  { theCollection[i].Print(); }
}



