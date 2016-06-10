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
// $Id: G4NormalNavigation.hh 90009 2015-05-08 07:42:39Z gcosmo $
//
// 
// class G4NormalNavigation
//
// Class description:
//
// Utility for navigation in volumes containing only G4PVPlacement
// daughter volumes.

// History:
// - Created. Paul Kent, Aug 96
// --------------------------------------------------------------------
#ifndef G4NORMALNAVIGATION_HH
#define G4NORMALNAVIGATION_HH

#include <iomanip>

#include "G4NavigationHistory.hh"

#include "G4VPhysicalVolume.hh"
#include "G4LogicalVolume.hh"
#include "G4VSolid.hh"
#include "G4ThreeVector.hh"
#include "G4AuxiliaryNavServices.hh"

class G4NavigationLogger;

class G4NormalNavigation
{
  public:  // with description

    G4NormalNavigation();
      // Constructor

    ~G4NormalNavigation();
      // Destructor

    inline G4bool LevelLocate( G4NavigationHistory &history,
                  const G4VPhysicalVolume *blockedVol,
                  const G4int blockedNum,
                  const G4ThreeVector &globalPoint,
                  const G4ThreeVector* globalDirection,
                  const G4bool pLocatedOnEdge, 
                        G4ThreeVector &localPoint);
      // Search positioned volumes in mother at current top level of history
      // for volume containing globalPoint. Do not test the blocked volume.
      // If a containing volume is found, `stack' the new volume and return
      // true, else return false (the point lying in the mother but not any
      // of the daughters). localPoint = global point in local system on entry,
      // point in new system on exit.

    G4double ComputeStep( const G4ThreeVector &localPoint,
                          const G4ThreeVector &localDirection,
                          const G4double currentProposedStepLength,
                                G4double &newSafety,
                                G4NavigationHistory &history,
                                G4bool &validExitNormal,
                                G4ThreeVector &exitNormal,
                                G4bool &exiting,
                                G4bool &entering,
                                G4VPhysicalVolume *(*pBlockedPhysical),
                                G4int &blockedReplicaNo );

    G4double ComputeSafety( const G4ThreeVector &globalpoint,
                            const G4NavigationHistory &history,
                            const G4double pMaxLength=DBL_MAX );

    G4int GetVerboseLevel() const;
    void  SetVerboseLevel(G4int level);
      // Get/Set Verbose(ness) level.
      // [if level>0 && G4VERBOSE, printout can occur]

    inline void  CheckMode(G4bool mode);
      // Run navigation in "check-mode", therefore using additional
      // verifications and more strict correctness conditions.
      // Is effective only with G4VERBOSE set.

  private:

    G4bool fCheck; 
    G4NavigationLogger* fLogger;
};

#include "G4NormalNavigation.icc"

#endif
