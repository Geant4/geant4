// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4SceneDataObjectList.hh,v 1.1 1999-01-07 16:15:17 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
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
