// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4VPersistentSubMan.hh,v 1.2 1999/11/25 11:33:07 morita Exp $
// GEANT4 tag $Name: geant4-01-00 $
//
// class G4VPersistentSubMan 
//
// Abstract submanager class to be used by G4PersistencyManager
//
// Member functions:
// =================
//  virtual void SetVerboseLevel(G4int verboseLevel)
//    Set verbose level
//  inline G4Pint GetVerboseLevel()
//    Get verbose level
//
// Member data:
// ============
//  G4int f_verboseLevel;
//    Internal flag for verbose message printing
//
// History:
// 98.10.30 Y.Morita  Splited from G4PersistencyManager

#ifndef G4VPersistentSubMan_hh
#define G4VPersistentSubMan_hh 1

#include "globals.hh"
#include "G4PersistentTypes.hh"

class G4VPersistentSubMan 
{
  friend class G4PersistencyManager;

  protected:
      // to be used by G4PersistencyManager only
      G4VPersistentSubMan();
      virtual ~G4VPersistentSubMan();

  protected:
      // interface with G4PersistencyManager
      inline virtual void SetVerboseLevel(G4int verboseLevel)
      { f_verboseLevel = verboseLevel; }

  public:
      inline G4Pint GetVerboseLevel()
      { return f_verboseLevel; };

  protected:
      G4int f_verboseLevel;

};

#endif

