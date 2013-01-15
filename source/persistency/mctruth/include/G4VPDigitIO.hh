//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
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
      static G4ThreadLocal G4VPDigitIO* f_G4VPDigitIO;
      G4DCIOcatalog*    f_catalog;

}; // End of class G4VPDigitIO

#endif

