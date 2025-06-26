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
// A class that provides an estimate of the isotropic safety, the minimum
// distance from a global point to the nearest boundary of the current volume
// or the nearest daughter volumes.
// This estimate can be an underestimate, either because a solid provides an
// underestimate (for speed) or in order to avoid substantial additional
// computations.
// Obtains from the navigator the current transformation history.

// Author: John Apostolakis (CERN), February 2023
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

/**
 * @brief G4SafetyCalculator is a class that provides an estimate of the
 * isotropic safety (the minimum distance from a global point to the nearest
 * boundary of the current volume or the nearest daughter volumes).
 */

class G4SafetyCalculator
{
  public:

    /**
     * Constructor, initialisers and setup.
     */
    G4SafetyCalculator( const G4Navigator& navigator,
                        const G4NavigationHistory& navHistory );   

    /**
     * Copy constructor & assignment operator not allowed.
     */
    G4SafetyCalculator(const G4SafetyCalculator&) = delete;
    G4SafetyCalculator& operator=(const G4SafetyCalculator&) = delete;

    /**
     * Destructor. No actions.
     */
    ~G4SafetyCalculator() = default;

    /**
     * Calculates the isotropic distance to the nearest boundary from the
     * specified point in the global coordinate system. 
     *  @param[in] globalPoint The point in global coordinates; it *must* be
     *             located exactly within the current volume (it also must
     *             *not* be in a daughter volume.
     *  @param[in] physicalVolume Current volume.
     *  @param[in] pProposedMaxLength The calculation will not look beyond
     *             the proposed maximum length to avoid extra volume safety
     *             calculations.
     *  @param[in] verbose Flag to enable verbosity (default is false).
     *  @returns An underestimate of the safety distance (and typically will be
     *           if complex volumes are involved.
     */
    G4double SafetyInCurrentVolume(const G4ThreeVector& globalpoint,
                                         G4VPhysicalVolume* physicalVolume,
                                   const G4double pProposedMaxLength = DBL_MAX,
                                         G4bool verbose = false );

    /**
     * Accessor & modifier for custom external navigation.
     */
    G4VExternalNavigation* GetExternalNavigation() const;
    void SetExternalNavigation(G4VExternalNavigation* externalNav);
   
    /**
     * Compares estimates of the safety, and reports if found difference(s).
     */
    void CompareSafetyValues( G4double oldSafety,
                              G4double newValue,
                              G4VPhysicalVolume* motherPhysical,
                        const G4ThreeVector &globalPoint,
                              G4bool keepState,
                              G4double maxLength,
                              G4bool enteredVolume,
                              G4bool exitedVolume );

  protected:

    /**
     * Prepare state of sub-navigators by informing them of current point.
     *  @param[in] pointLocal Point in local coordinates.
     *  @param[in] motherPhysical Pointer to current volume where to relocate
     *             in case an external custom navigator is used.
     */
    void QuickLocateWithinVolume(const G4ThreeVector& pointLocal,
                                       G4VPhysicalVolume* motherPhysical);
   
    /**
     * Computes point in local coordinates system, given a position
     * vector in world coordinate system.
     *  @param[in] rGlobPoint Point in global coordinates.
     *  @returns The point in local coordinates system.
     */
    inline G4ThreeVector ComputeLocalPoint(const G4ThreeVector& rGlobPoint) const;

    /**
     * Computes the local direction of the specified vector in the reference
     * system of the volume that was found by LocateGlobalPointAndSetup().
     *  @param[in] pVec Vector in global coordinates.
     *  @returns The local direction of the specified vector in global
     *           coordinates.
     */
    inline G4ThreeVector ComputeLocalAxis(const G4ThreeVector& pVec) const;

    /**
     * Characterises the daughter of logical volume.
     */
    inline EVolume CharacteriseDaughters(const G4LogicalVolume* pLog) const;

    /**
     * Gets regular structure ID of first daughter.
     */
    inline G4int GetDaughtersRegularStructureId(const G4LogicalVolume* pLv) const;

  private:

    // BEGIN -- Tracking Invariants part 1 ------------------------------------

    /** Associated navigator. Needed for optimisation. */
    const G4Navigator& fNavigator;

    /** Associated navigator's navigation history. Transformation and history
        of the current path through the geometrical hierarchy. */
    const G4NavigationHistory& fNavHistory;

    // END   -- Tracking Invariants part 1 ------------------------------------

    /** Cached tolerance. */
    G4double fkCarTolerance; 
   
    // BEGIN State information ------------------------------------------------

    /** Previous safety origin. */
    G4ThreeVector fPreviousSftOrigin;

    /** Memory of last safety origin & value. Used in ComputeStep() to ensure
        that origin of current Step is in the same volume as the point of the
        last relocation. */
    G4double fPreviousSafety = 0.0; 

    /** Helpers/Utility classes - their state can change. */
    G4NormalNavigation fnormalNav;
    G4VoxelNavigation fvoxelNav;
    G4ParameterisedNavigation fparamNav;
    G4ReplicaNavigation freplicaNav;
    G4RegularNavigation fregularNav;
    G4VExternalNavigation* fpExternalNav = nullptr;
    G4VoxelSafety fVoxelSafety;
};

// Auxiliary inline methods -- copied from G4Navigator
//
#include "G4SafetyCalculator.icc"

#endif
