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
// $Id: G4PersistentEventMan.hh,v 1.13 2001/07/11 10:02:15 gunter Exp $
// GEANT4 tag $Name: geant4-04-00 $
//

// Class Description:
//   A Utility class to be used by G4PersistencyManager.
// Average users do not need to use this class.
//
// This class is not persistent-capable. 
//

// Member functions:
// =================
//  G4bool Store(HepDbApplication* dbApp,
//               const G4Event* anEvent);
//    Store anEvent and associated objects into database.
//  G4bool Retrieve(HepDbApplication* dbApp,
//                  G4Event*& anEvent);
//    Retrieve anEvent and associated objects from database.
//
// History:
// 98.10.30 Y.Morita  Splited from G4PersistencyManager

#ifndef G4PERSISTENTEVENTMAN_HH
#define G4PERSISTENTEVENTMAN_HH 1

#include "globals.hh"
#include "HepODBMS/odbms/HepODBMS.h"

#include "G4VPersistentSubMan.hh"
#include "G4VPersistentSubDbMan.hh"

#include "G4PEvent.hh"

class HepDbApplication;
class G4Event;
class G4PersistentHitMan;
class G4PersistentDigitMan;

class G4PersistentEventMan 
 : public G4VPersistentSubMan, public G4VPersistentSubDbMan
{
  friend class G4PersistencyManager;

  private:
      // to be used by G4PersistencyManager only
      G4PersistentEventMan();
      G4PersistentEventMan( G4PersistentHitMan*   aPHCMan,
                            G4PersistentDigitMan* aPDCMan );
      ~G4PersistentEventMan();

  private:
      G4PersistentHitMan*   f_PHCMan;
      G4PersistentDigitMan* f_PDCMan;
      G4int         f_currentEventID;
      HepRef(G4PEvent) f_currentPEvent;

  public: // With description
      inline G4int CurrentEventID()
      { return f_currentEventID; }
        // Returns the current event id

  private:
      G4bool Store( HepDbApplication* dbApp,
                    const G4Event* anEvent );
      G4bool Retrieve( HepDbApplication* dbApp,
                       G4Event*& anEvent );

};

#endif

