// $Id: G4VPDigitIO.hh,v 1.3 2002-12-04 13:57:29 morita Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// File: G4VPDigitIO.hh
//
// History:
//   '01.08.10  Youhei Morita  Initial creation (with "fadsclass3")

#ifndef V_P_DIGIT_I_O_HH
#define V_P_DIGIT_I_O_HH 1

#include "G4DCofThisEvent.hh"
#include "G4DCIOcatalog.hh"
#include "G4VPDigitsCollectionIO.hh"

// Class Description:
//   Abstract base class for storing and retrieving digit collections

class G4VPDigitIO
{
    public: // With description
      G4VPDigitIO();
      // Constructor

      virtual ~G4VPDigitIO() {};
      // Destructor

    public: // With description
      G4VPDigitIO* GetVPDigitIO() { return f_G4VPDigitIO; };
      // Returns the pointer of the digit collection I/O manager.

      virtual G4bool Store(const G4DCofThisEvent*) =0;
      // Pure virtual method for storing digit collections of this event.
      // Each persistency package should implement a concrete method
      // of storing the digit collection of this event with this signature.

      virtual G4bool Retrieve(G4DCofThisEvent*&) =0;
      // Pure virtual method for retrieving digit collections of this event.
      // Each persistency package should implement a concrete method
      // of storing the digit collection of this event with this signature.

      void SetVerboseLevel(int v);
      // Set verbose level.

    protected:
      void SetG4VPDigitIO(G4VPDigitIO* digitMan) { f_G4VPDigitIO = digitMan; };
      // Registers the digit collection I/O manager.

    protected:
      G4int m_verbose;
      static G4VPDigitIO* f_G4VPDigitIO;
      G4DCIOcatalog*    f_catalog;

}; // End of class G4VPDigitIO

#endif

