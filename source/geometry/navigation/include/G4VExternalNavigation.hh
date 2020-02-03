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
// G4VExternalNavigation
//
// Class description:
//
// Pure virtual class to be specialised by the user for tracking with
// an external navigation

// Authors: V.Vlachoudis, G.Cosmo - CERN, 2019
// --------------------------------------------------------------------
#ifndef G4VEXTERNALNAVIGATION_HH
#define G4VEXTERNALNAVIGATION_HH

#include "G4NavigationHistory.hh"
#include "G4VPhysicalVolume.hh"
#include "G4LogicalVolume.hh"
#include "G4VSolid.hh"
#include "G4ThreeVector.hh"

class G4VExternalNavigation
{
  public:  // with description

    G4VExternalNavigation();
      // Constructor
   
    virtual ~G4VExternalNavigation();
      // Destructor
   
    virtual G4bool LevelLocate( G4NavigationHistory& history,
                                const G4VPhysicalVolume* blockedVol,
                                const G4int blockedNum,
                                const G4ThreeVector& globalPoint,
                                const G4ThreeVector* globalDirection,
                                const G4bool pLocatedOnEdge,
                                G4ThreeVector& localPoint) = 0;
      // Search positioned volumes in mother at current top level of history
      // for volume containing globalPoint. Do not test the blocked volume.
      // If a containing volume is found, `stack' the new volume and return
      // true, else return false (the point lying in the mother but not any
      // of the daughters). localPoint = global point in local system on entry,
      // point in new system on exit.

    virtual G4double ComputeStep( const G4ThreeVector& localPoint,
                                  const G4ThreeVector& localDirection,
                                  const G4double currentProposedStepLength,
                                  G4double& newSafety,
                                  G4NavigationHistory& history,
                                  G4bool& validExitNormal,
                                  G4ThreeVector& exitNormal,
                                  G4bool& exiting,
                                  G4bool& entering,
                                  G4VPhysicalVolume** pBlockedPhysical,
                                  G4int& blockedReplicaNo ) = 0;
     // Compute the length of a step to the next boundary.
     // Ignore the (input) pBlockedPhysical volume (with replica/parameterisation
     //   number 'blockedReplica')
     // Identify the next candidate volume (if a daughter of current volume),
     //   and return it in pBlockedPhysical, blockedReplicaNo
     // In/Out  Navigation history: to be update for volume of next intersection.
     // In/out 'newsafety' is the known isotropic safety of the initial point:
     //     an earlier (likely crude) estimate on input,
     //     to be updated with new estimate for the same (start) point.
   
    virtual G4double ComputeSafety( const G4ThreeVector& globalpoint,
                                    const G4NavigationHistory& history,
                                    const G4double pMaxLength = DBL_MAX ) = 0;

    virtual G4VExternalNavigation* Clone() = 0;

    // Optional methods - may be necessary under particular circumstances

    virtual EInside Inside( const G4VSolid*      solid,
                            const G4ThreeVector& position,
                            const G4ThreeVector& direction );
     // Special 'Inside' call that includes direction of next motion
     //   provided for potential optimisations.

    virtual void RelocateWithinVolume( G4VPhysicalVolume*  motherPhysical,
                                       const G4ThreeVector& localPoint );
     //   Update any relevant internal state to take account that
     //      - the location has been moved to 'localPoint'
     //      - it remains in the current (mother) physical volume 'motherPhysical'
   
    inline G4int GetVerboseLevel() const { return fVerbose; }
    inline void  SetVerboseLevel(G4int level) { fVerbose = level; }
      // Get/Set verbosity level.

    inline void  CheckMode(G4bool mode) { fCheck = mode; }
      // Run navigation in "check-mode", therefore using additional
      // verifications and more strict correctness conditions.
      // Should be effective only with G4VERBOSE set.

  protected:

    G4bool fCheck = false;
    G4int fVerbose = 0;
};

#endif
