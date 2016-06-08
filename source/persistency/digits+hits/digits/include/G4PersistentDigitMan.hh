// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4PersistentDigitMan.hh,v 1.7 2000/12/15 08:04:13 morita Exp $
// GEANT4 tag $Name: geant4-03-01 $
//

// Class Description:
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
  public: // with description
      static G4PersistentDigitMan* GetPersistentDigitMan();
      // returns a pointer to singleton instance of this class.
      static G4PersistentDigitMan* get_PersistentDigitManIfExist();
      // returns a pointer to singleton instance of this class.

  public:
      G4PersistentDigitMan();
      ~G4PersistentDigitMan();

  private: 
      static G4PersistentDigitMan * f_PersistentDigitMan;
      HepRef(G4PDCofThisEvent) f_CurrentPDCofThisEvent;

  public: // with description
      inline void SetCurrentPDCofThisEvent( HepRef(G4PDCofThisEvent) aPDC)
      { f_CurrentPDCofThisEvent = aPDC; };
      // sets a smart pointer of the collection of digits collections
      // of this event.  This method should be invoked by user
      // sensitive detector only once per event.
      inline HepRef(G4PDCofThisEvent) GetCurrentPDCofThisEvent()
      { return f_CurrentPDCofThisEvent; };
      // returns a smart pointer of the collection of digits collections
      // of this event.

};

#endif

