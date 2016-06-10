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
// $Id: G4GlobalMagFieldMessenger.hh 66536 2012-12-19 14:32:36Z ihrivnac $
//
// class G4GlobalMagFieldMessenger
//
// Class description:
//
// Global uniform magnetic field messenger class. 
//
// It defines UI commands:
// - /globalField/setValue vx vy vz unit
// - /globalField/verbose level
//
// It creates/deletes the global uniform magnetic field and
// activates/inactivates it according to the set field value.
// The field value can be changed either interactively via 
// the UI command or via SetFieldValue() function.

// Author: Ivana Hrivnacova, 28/08/2013  (ivana@ipno.in2p3.fr)
// --------------------------------------------------------------------
#ifndef G4GlobalMagFieldMessenger_hh
#define G4GlobalMagFieldMessenger_hh 1

#include "globals.hh"
#include "G4UImessenger.hh"

class G4UniformMagField;
class G4UIdirectory;
class G4UIcmdWith3VectorAndUnit;
class G4UIcmdWithAnInteger;

class G4GlobalMagFieldMessenger : public G4UImessenger
{
  public:  // with description

    G4GlobalMagFieldMessenger(const G4ThreeVector& value = G4ThreeVector());
    virtual ~G4GlobalMagFieldMessenger();
    
    virtual void SetNewValue(G4UIcommand*, G4String);

    void  SetFieldValue(const G4ThreeVector& value);
    G4ThreeVector GetFieldValue() const;
    
    inline void  SetVerboseLevel(G4int verboseLevel);
    inline G4int GetVerboseLevel() const;
    
  private:

    void SetField(const G4ThreeVector& value, const G4String& inFunction);
    
    G4UniformMagField*  fMagField;
    G4int               fVerboseLevel;

    G4UIdirectory*      fDirectory;
    G4UIcmdWith3VectorAndUnit* fSetValueCmd;
    G4UIcmdWithAnInteger*      fSetVerboseCmd;
};

// inline functions

inline void  G4GlobalMagFieldMessenger::SetVerboseLevel(G4int verboseLevel)
{ fVerboseLevel = verboseLevel; }

inline G4int G4GlobalMagFieldMessenger::GetVerboseLevel() const
{ return fVerboseLevel; }
    
#endif
