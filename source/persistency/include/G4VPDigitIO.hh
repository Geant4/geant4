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

