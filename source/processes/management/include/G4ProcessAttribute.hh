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
// G4ProcessAttribute
//
// Class description:
//
// This class is used exclusively by G4ProcessManager for booking.

// Author: H.Kurashige, 2 December 1997
// --------------------------------------------------------------------
#ifndef G4ProcessAttribute_hh
#define G4ProcessAttribute_hh 1

#include "globals.hh"
#include "G4ios.hh"

#include "G4ProcessManager.hh"

class G4VProcess;

class G4ProcessAttribute
{
  friend class G4ProcessManager;

  public:

    G4ProcessAttribute();
    G4ProcessAttribute(const G4VProcess* aProcess);
    G4ProcessAttribute(const G4ProcessAttribute& right);
      // Constructors

    ~G4ProcessAttribute();
      // Destructor

    G4ProcessAttribute& operator=(const G4ProcessAttribute& right);
      // Assignment operator

    inline G4bool operator==(const G4ProcessAttribute &right) const;
    inline G4bool operator!=(const G4ProcessAttribute &right) const;
      // Equality operators

  protected:

    G4VProcess* pProcess = nullptr;
      // Pointer to G4VProcess

    G4bool isActive = true;
      // Flag for activation/inactivation

    G4int idxProcessList = -1;
      // Index to a ProcessVector for theProcessList

    G4int idxProcVector[G4ProcessManager::SizeOfProcVectorArray];
      // Index to ProcessVectors for Doit() and GetPhysicalInteractionLength()
      // methods. A value of -1 means "not applicable"

    G4int ordProcVector[G4ProcessManager::SizeOfProcVectorArray];
      // Ordering parameter 
};

// ------------------------
// Inline methods
// ------------------------

inline 
G4bool G4ProcessAttribute::operator==(const G4ProcessAttribute& right) const
{
  return this->pProcess == right.pProcess;
}

inline 
G4bool G4ProcessAttribute::operator!=(const G4ProcessAttribute& right) const
{
  return this->pProcess != right.pProcess;
}

#endif
