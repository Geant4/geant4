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
// $Id: G4PersistentHitMan.hh,v 1.6.4.1 2001/06/28 19:11:25 gunter Exp $
// GEANT4 tag $Name:  $
//

// Class Description:
//
// A Class responsible for registring HitsCollection
// during event tracking.
//
// The class is `singleton', with access via
//   G4PersistentHitMan::GetPersistentHitMan.
//
// This class itself is not persistent-capable. 
//

// Member functions:
// =================
//  static G4PersistentHitMan* GetPersistentHitMan();
//  static G4PersistentHitMan* get_PersistentHitManIfExist();
//   Return ptr to singleton instance of the class.
//
//  G4PersistentHitMan();
//  ~G4PersistentHitMan();
//
// Member data:
// ============
// static G4PersistentHitMan* fPersistentHitMan
//   Ptr to the unique instance of class
//
// History:
// 99.11.16 Y.Morita  Initial version

#ifndef G4PersistentHitMan_DDL
#define G4PersistentHitMan_DDL 1

#include "HepODBMS/odbms/HepODBMS.h"

#include "G4VPersistentSubMan.hh"
#include "G4VPersistentSubDbMan.hh"

#include "G4PHCofThisEvent.hh"

class G4PersistentHitMan 
 : public G4VPersistentSubMan, public G4VPersistentSubDbMan
{
  public: // with description
      static G4PersistentHitMan* GetPersistentHitMan();
      // returns a pointer to singleton instance of this class.
      static G4PersistentHitMan* get_PersistentHitManIfExist();
      // returns a pointer to singleton instance of this class.

  public:
      G4PersistentHitMan();
      ~G4PersistentHitMan();

  private: 
      static G4PersistentHitMan * f_PersistentHitMan;
      HepRef(G4PHCofThisEvent) f_CurrentPHCofThisEvent;

  public: // with description
      inline void SetCurrentPHCofThisEvent( HepRef(G4PHCofThisEvent) aPHC)
      { f_CurrentPHCofThisEvent = aPHC; };
      // sets a smart pointer of the collection of hits collections
      // of this event.  This method should be invoked by user
      // sensitive detector only once per event.
      inline HepRef(G4PHCofThisEvent) GetCurrentPHCofThisEvent()
      { return f_CurrentPHCofThisEvent; };
      // returns a smart pointer of the collection of hits collections
      // of this event.

};

#endif

