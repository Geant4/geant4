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
// G4SafetyCalculator
//
// Class description:
//
// A class that provides an estimate of the isotropic safety - the
//   minimum distance from a global point to the nearest boundary
//   of the current volume or the nearest daughter volumes.
// This estimate can be an underestimate, either because a solid
//   provides an underestimate (for speed) or in order to avoid
//   substantial additional computations.
//
// Obtains from the navigator the current transformation history.

// Author: John Apostolakis, CERN - February 2023
// --------------------------------------------------------------------
#ifndef G4SafetyCalculator_HH
#define G4SafetyCalculator_HH 1

#include "geomdefs.hh"

#include "G4ThreeVector.hh"
#include "G4AffineTransform.hh"
#include "G4RotationMatrix.hh"

#include "G4LogicalVolume.hh"             // Used in inline methods
#include "G4TouchableHistoryHandle.hh"

#include "G4NavigationHistory.hh"
#include "G4NormalNavigation.hh"
#include "G4VoxelNavigation.hh"
#include "G4ParameterisedNavigation.hh"
#include "G4ReplicaNavigation.hh"
#include "G4RegularNavigation.hh"
#include "G4VExternalNavigation.hh"

#include "G4VoxelSafety.hh"

#include <iostream>

class G4VPhysicalVolume;

class G4SafetyCalculator
{
  public:

    G4SafetyCalculator( const G4Navigator& navigator,
                        const G4NavigationHistory& navHistory );   
      // Constructor - initialisers and setup.

    G4SafetyCalculator(const G4SafetyCalculator&) = delete;
    G4SafetyCalculator& operator=(const G4SafetyCalculator&) = delete;
      // Copy constructor & assignment operator not allowed.

    ~G4SafetyCalculator() = default;
      // Destructor. No actions.

    G4double SafetyInCurrentVolume(const G4ThreeVector& globalpoint,
                                         G4VPhysicalVolume* physicalVolume,
                                   const G4double pProposedMaxLength = DBL_MAX,
                                         G4bool verbose = false );
      // Calculate the isotropic distance to the nearest boundary from the
      // specified point in the global coordinate system. 
      // The globalpoint utilised *must* be located exactly within the
      // current volume (it also must *not* be in a daughter volume).
      // The value returned can be an underestimate (and typically will be
      // if complex volumes are involved).  
      // The calculation will not look beyond the proposed maximum length
      // to avoid extra volume safety calculations. The geometry must be closed.

    G4VExternalNavigation* GetExternalNavigation() const;
    void SetExternalNavigation(G4VExternalNavigation* externalNav);
      // Accessor & modifier for custom external navigation.
   
    void CompareSafetyValues( G4double oldSafety,
                              G4double newValue,
                              G4VPhysicalVolume* motherPhysical,
                        const G4ThreeVector &globalPoint,
                              G4bool keepState,
                              G4double maxLength,
                              G4bool enteredVolume,
                              G4bool exitedVolume );
      // Compare estimates of the safety, and report if difference(s) found.

  protected:

    void QuickLocateWithinVolume(const G4ThreeVector& pointLocal,
                                       G4VPhysicalVolume*  motherPhysical);
      // Prepare state of sub-navigators by informing them of current point.
   
    inline G4ThreeVector ComputeLocalPoint(const G4ThreeVector& rGlobPoint) const;
      // Return position vector in local coordinate system, given a position
      // vector in world coordinate system.

    inline G4ThreeVector ComputeLocalAxis(const G4ThreeVector& pVec) const;
      // Return the local direction of the specified vector in the reference
      // system of the volume that was found by LocalGlobalPointAndSetup().
      // The Local Coordinates of point in world coordinate system.

    inline EVolume CharacteriseDaughters(const G4LogicalVolume* pLog) const;
      // Characterise daughter of logical volume.

    inline G4int GetDaughtersRegularStructureId(const G4LogicalVolume* pLv) const;
      // Get regular structure ID of first daughter.

  private:

    // BEGIN -- Tracking Invariants part 1
    //
    const G4Navigator& fNavigator;
      // Associated navigator. Needed for details of current state,
      // for optimisation

    const G4NavigationHistory& fNavHistory;
      // Associated navigator's navigation history. Transformation and history
      // of the current path through the geometrical hierarchy.
    //
    // END   -- Tracking Invariants part 1

    G4double fkCarTolerance; 
      // Cached tolerance.
   
    // BEGIN State information
    //
    G4ThreeVector  fPreviousSftOrigin;
    G4double       fPreviousSafety = 0.0; 
      // Memory of last safety origin & value. Used in ComputeStep to ensure
      // that origin of current Step is in the same volume as the point of the
      // last relocation.

    // Helpers/Utility classes - their state can change
    //
    G4NormalNavigation fnormalNav;
    G4VoxelNavigation fvoxelNav;
    G4ParameterisedNavigation fparamNav;
    G4ReplicaNavigation freplicaNav;
    G4RegularNavigation fregularNav;
    G4VExternalNavigation* fpExternalNav = nullptr;
    G4VoxelSafety fVoxelSafety;
};

// Auxiliary inline methods -- copied from G4Navigator

// Return  local coordinates given point in the world coord system.
//
inline G4ThreeVector
G4SafetyCalculator::ComputeLocalPoint(const G4ThreeVector& pGlobalPoint) const
{
  return fNavHistory.GetTopTransform().TransformPoint(pGlobalPoint);
}

// Returns local direction given vector direction in world coord system.
//
inline G4ThreeVector
G4SafetyCalculator::ComputeLocalAxis(const G4ThreeVector& pVec) const
{
  return fNavHistory.GetTopTransform().TransformAxis(pVec);
}

inline EVolume
G4SafetyCalculator::CharacteriseDaughters(const G4LogicalVolume* pLog) const
{
  return pLog->CharacteriseDaughters();
}

inline G4int
G4SafetyCalculator::GetDaughtersRegularStructureId(const G4LogicalVolume* pLog) const
{
  G4int regId = 0;
  G4VPhysicalVolume *pVol;

  if ( pLog->GetNoDaughters() == 1 )
  {
    pVol = pLog->GetDaughter(0);
    regId = pVol->GetRegularStructureId();
  }
  return regId;
}

#endif
