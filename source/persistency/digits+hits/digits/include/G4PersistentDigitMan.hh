// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4PersistentDigitMan.hh,v 1.6 1999/11/29 20:11:16 morita Exp $
// GEANT4 tag $Name: geant4-01-00 $
//
// class G4PersistentDigitMan 
//
// A Class responsible for registring DigitsCollection
// during event tracking.
//
// The class is `singleton', with access via
//   G4PersistentDigitMan::GetPersistentDigitMan.
//
// This class itself is not persistent-capable. 
//
// Member functions:
// =================
//  static G4PersistentDigitMan* GetPersistentDigitMan();
//  static G4PersistentDigitMan* get_PersistentDigitManIfExist();
//   Return ptr to singleton instance of the class.
//
//  G4PersistentDigitMan();
//  ~G4PersistentDigitMan();
//
// Member data:
// ============
// static G4PersistentDigitMan* fPersistentDigitMan
//   Ptr to the unique instance of class
//
// History:
// 99.11.22 Y.Morita  Initial version

#ifndef G4PersistentDigitMan_DDL
#define G4PersistentDigitMan_DDL 1

#include "HepODBMS/odbms/HepODBMS.h"

#include "G4VPersistentSubMan.hh"
#include "G4VPersistentSubDbMan.hh"

#include "G4PDCofThisEvent.hh"

class G4PersistentDigitMan 
 : public G4VPersistentSubMan, public G4VPersistentSubDbMan
{
  public:
      static G4PersistentDigitMan* GetPersistentDigitMan();
      static G4PersistentDigitMan* get_PersistentDigitManIfExist();

      G4PersistentDigitMan();
      ~G4PersistentDigitMan();

  private: 
      static G4PersistentDigitMan * f_PersistentDigitMan;
      HepRef(G4PDCofThisEvent) f_CurrentPDCofThisEvent;

  public:
      inline void SetCurrentPDCofThisEvent( HepRef(G4PDCofThisEvent) aPDC)
      { f_CurrentPDCofThisEvent = aPDC; };
      inline HepRef(G4PDCofThisEvent) GetCurrentPDCofThisEvent()
      { return f_CurrentPDCofThisEvent; };

};

#endif

