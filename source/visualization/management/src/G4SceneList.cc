// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4SceneList.cc,v 1.1 1999-01-09 16:31:14 allison Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// John Allison  9th August 1998

#include "G4SceneList.hh"

G4SceneList::G4SceneList
(
 unsigned (*hashFun)(const G4String&),      // Hashing function
 size_t size                                // No. of buckets
 ):
RWTValHashDictionary <G4String, G4Scene> (hashFun, size)
{}

G4SceneListIterator::G4SceneListIterator
(G4SceneList& list):
RWTValHashDictionaryIterator <G4String, G4Scene> (list)
{}
