// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: PersEx02TrackerHitsCollection.cc,v 1.3 1999/11/29 18:33:29 morita Exp $
// GEANT4 tag $Name: geant4-01-00 $
//

#include "PersEx02TrackerHitsCollection.hh"


PersEx02TrackerHitsCollection::PersEx02TrackerHitsCollection
                                    (G4String dName, G4String aName)
 : G4PVHitsCollection(dName, aName), elems(0)
{;}

PersEx02TrackerHitsCollection::~PersEx02TrackerHitsCollection()
{;}

void PersEx02TrackerHitsCollection::DrawAllHits()
{;}

void PersEx02TrackerHitsCollection::PrintAllHits()
{;}

