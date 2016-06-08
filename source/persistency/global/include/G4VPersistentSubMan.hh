//
// ********************************************************************
// * DISCLAIMER                                                       *
// *                                                                  *
// * The following disclaimer summarizes all the specific disclaimers *
// * of contributors to this software. The specific disclaimers,which *
// * govern, are listed with their locations in:                      *
// *   http://cern.ch/geant4/license                                  *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.                                                             *
// *                                                                  *
// * This  code  implementation is the  intellectual property  of the *
// * GEANT4 collaboration.                                            *
// * By copying,  distributing  or modifying the Program (or any work *
// * based  on  the Program)  you indicate  your  acceptance of  this *
// * statement, and all its terms.                                    *
// ********************************************************************
//
//
// $Id: G4VPersistentSubMan.hh,v 1.5 2001/07/11 10:02:25 gunter Exp $
// GEANT4 tag $Name: geant4-04-01 $

// Class Description:
//   Abstract submanager class to be used by G4PersistencyManager

// Member functions:
// =================
//  virtual void SetVerboseLevel(G4int verboseLevel)
//    Set verbose level
//  inline G4int GetVerboseLevel()
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
        // Set the verbose level

  public: // With description
      inline G4int GetVerboseLevel()
      { return f_verboseLevel; };
        // Get the current verbose level

  protected:
      G4int f_verboseLevel;

};

#endif

