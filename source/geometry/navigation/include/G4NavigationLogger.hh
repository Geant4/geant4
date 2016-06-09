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
//
// $Id$
//
// 
// class G4NavigationLogger
//
// Class description:
//
// Simple utility class for use by navigation systems
// for verbosity and check-mode.

// History:
// - Created. Gabriele Cosmo, November 2010
// --------------------------------------------------------------------
#ifndef G4NAVIGATIONLOGGER_HH
#define G4NAVIGATIONLOGGER_HH

#include "G4NavigationHistory.hh"
#include "G4VPhysicalVolume.hh"
#include "G4LogicalVolume.hh"
#include "G4VSolid.hh"
#include "G4ThreeVector.hh"

class G4NavigationLogger
{
  public:  // with description

    G4NavigationLogger(const G4String& id);
   ~G4NavigationLogger();

    void PreComputeStepLog  (const G4VPhysicalVolume* motherPhysical,
                                   G4double motherSafety,
                             const G4ThreeVector& localPoint) const;
    void AlongComputeStepLog(const G4VSolid* sampleSolid,
                             const G4ThreeVector& samplePoint,
                             const G4ThreeVector& sampleDirection,
                             const G4ThreeVector& localDirection,
                                   G4double sampleSafety,
                                   G4double sampleStep) const;
    void PostComputeStepLog (const G4VSolid* motherSolid,
                             const G4ThreeVector& localPoint,
                             const G4ThreeVector& localDirection,
                                   G4double motherStep,
                                   G4double motherSafety) const;
    void ComputeSafetyLog   (const G4VSolid* solid,
                             const G4ThreeVector& point,
                                   G4double safety,
                                   G4bool banner) const;
    void PrintDaughterLog   (const G4VSolid* sampleSolid,
                             const G4ThreeVector& samplePoint,
                                   G4double sampleSafety,
                                   G4double sampleStep) const;
  public:  // without description

    inline G4int GetVerboseLevel() const  { return fVerbose; }
    inline void  SetVerboseLevel(G4int level)  { fVerbose = level; }

  private:

    G4String fId;     // Navigation type
    G4int fVerbose;   // Verbosity level
};

#endif
