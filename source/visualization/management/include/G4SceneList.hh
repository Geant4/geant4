// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4SceneList.hh,v 1.2 1999-01-09 16:30:43 allison Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// John Allison  9th August 1998

#ifndef G4SCENELIST_HH
#define G4SCENELIST_HH

#include "G4Scene.hh"
#include "globals.hh"
#include <rw/tvhdict.h>

class G4SceneList:
  public RWTValHashDictionary <G4String, G4Scene>
// Each scene is keyed to a name.
{
public:
  G4SceneList
  (
   unsigned (*)(const G4String&),      // Hashing function
   size_t size = RWDEFAULT_CAPACITY    // No. of buckets
   );
  friend class G4SceneListIterator;
};

class G4SceneListIterator:
  public RWTValHashDictionaryIterator <G4String, G4Scene>
{
public:
  G4SceneListIterator (G4SceneList&);
};

#endif
