// $Id: G4DigitRootIO.hh,v 1.2 2002-12-04 14:12:26 morita Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// File: G4DigitRootIO.hh
//
// History:
//   '02.5.7  Youhei Morita  Initial creation

#ifndef DIGIT_ROOT_IO_HH
#define DIGIT_ROOT_IO_HH 1

#include "G4DCofThisEvent.hh"
#include "G4PersistencyCenter.hh"
#include "G4RootIOManager.hh"
#include "G4DCIOcatalog.hh"

// Class inherited:
#include "G4VPDigitIO.hh"

// Class Description:
//   Manager class to store and retrieve Digit objects.
// 
//   This is a singleton class and should be constructed only
//   by GetDigitRootIO().

class G4DigitRootIO
 : public G4VPDigitIO
{
    public: // With description
      G4DigitRootIO();
      // Constructor

      virtual ~G4DigitRootIO();
      // Destructor

    public: // With description
      static G4DigitRootIO* GetDigitRootIO();
      // Construct a new singleton G4DigitRootIO object if it does not exist.

      bool Store(const G4DCofThisEvent* dcevt);
      // Store digit collections.
      // Concrete class of RootVDigitsCollectionIO must be registered
      // with G4VPDigitIO::AddDCIOmanager() before calling this method.

      bool Retrieve(G4DCofThisEvent*& dcevt);
      // Retrieve digit collections
      // Concrete class of RootVDigitsCollectionIO must be registered
      // with G4VPDigitIO::AddDCIOmanager() before calling this method.
      // Also a pointer of persistent digit collection must be set with
      // SetCurrentPDCofThisEvent().

}; // End of class G4DigitRootIO

#endif

