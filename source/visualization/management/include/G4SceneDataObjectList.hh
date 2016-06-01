// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4SceneDataObjectList.hh,v 2.1 1998/08/09 17:54:15 allison Exp $
// GEANT4 tag $Name: geant4-00 $
//
// 
// John Allison  9th August 1998

#ifndef G4SCENEDATAOBJECTLIST_HH
#define G4SCENEDATAOBJECTLIST_HH

#include "G4SceneData.hh"
#include "globals.hh"
#include <rw/tvhdict.h>

class G4SceneDataObjectList:
  public RWTValHashDictionary <G4String, G4SceneData>
// Each scene data object is keyed to a name.
{
public:
  G4SceneDataObjectList
  (
   unsigned (*)(const G4String&),      // Hashing function
   size_t size = RWDEFAULT_CAPACITY    // No. of buckets
   );
  friend class G4SceneDataObjectListIterator;
};

class G4SceneDataObjectListIterator:
  public RWTValHashDictionaryIterator <G4String, G4SceneData>
{
public:
  G4SceneDataObjectListIterator (G4SceneDataObjectList&);
};

#endif
