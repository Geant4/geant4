// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4SceneDataObjectList.cc,v 2.2 1998/08/10 10:19:09 allison Exp $
// GEANT4 tag $Name: geant4-00 $
//
// 
// John Allison  9th August 1998

#include "G4SceneDataObjectList.hh"

G4SceneDataObjectList::G4SceneDataObjectList
(
 unsigned (*hashFun)(const G4String&),      // Hashing function
 size_t size                                // No. of buckets
 ):
RWTValHashDictionary <G4String, G4SceneData> (hashFun, size)
{}

G4SceneDataObjectListIterator::G4SceneDataObjectListIterator
(G4SceneDataObjectList& list):
RWTValHashDictionaryIterator <G4String, G4SceneData> (list)
{}
