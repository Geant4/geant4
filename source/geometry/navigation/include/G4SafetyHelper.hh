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
// G4SafetyHelper
//
// Class description:
//
// This class is a helper for physics processes which require 
// knowledge of the safety, and the step size for the 'mass' geometry.

// Author: John Apostolakis (CERN), 5 July 2006
// --------------------------------------------------------------------
#ifndef G4SAFETYHELPER_HH
#define G4SAFETYHELPER_HH 1

#include <vector>

#include "G4Types.hh"
#include "G4ThreeVector.hh"
#include "G4Navigator.hh"

class G4PathFinder;

/**
 * @brief G4SafetyHelper is a helper class for physics processes which require 
 * knowledge of the safety, and the step size for the 'mass' geometry.
 */

class G4SafetyHelper
{
  public:

    /**
     * Constructor and default Destructor.
     */
    G4SafetyHelper(); 
   ~G4SafetyHelper() = default;

    /**
     * Computes the distance in the mass geometry.
     *  @param[in] position Point in global coordinates.
     *  @param[in] direction Direction.
     *  @param[in] currentMaxStep Proposed step length to nearest boundary.
     *  @param[in,out] newSafety New safety.
     *  @returns The linear step for mass geometry.
     */
    G4double CheckNextStep( const G4ThreeVector& position, 
                            const G4ThreeVector& direction,
                            const G4double currentMaxStep,
                                  G4double& newSafety );

    /**
     * Computes the safety distance for all geometries.
     *  @param[in] pGlobalPoint Point in global coordinates.
     *  @param[in] maxRadius Radius of interest (e.g. maximum displacement).
     *             Giving this, one can reduce the average computational
     *             cost. If not provided, the real isotropic safety is computed.
     *  @returns The safety distance for all geometries.
     */
    G4double ComputeSafety( const G4ThreeVector& pGlobalPoint,
                            G4double maxRadius = DBL_MAX );

    /**
     * Locates the point for all geometries.
     *  @param[in] pGlobalPoint Point in global coordinates.
     *  @param[in] direction Direction.
     */
    void Locate(const G4ThreeVector& pGlobalPoint,
                const G4ThreeVector& direction);

    /**
     * Relocates the point in the volume of interest.
     *  @param[in] pGlobalPoint Point in global coordinates.
     */
    void ReLocateWithinVolume(const G4ThreeVector& pGlobalPoint );

    /**
     * Enables navigation in parallel geometries.
     *  @param[in] parallel Flag to have parallel worlds considered.
     *             Alternative is to use single (mass) navigator directly.
     */
    inline void EnableParallelNavigation(G4bool parallel);

    /**
     * Checks for new navigator for tracking, and reinitialises pointer.
     */
    void InitialiseNavigator();

    /**
     * Verbosity control.
     *  @param[in] lev The new verbosity level to enable.
     *  @returns The old verbosity level.
     */
    inline G4int SetVerboseLevel( G4int lev );

    /**
     * Retrieves the world volume of the mass geometry.
     *  @returns The pointer to the mass geometry world volume.
     */
    inline G4VPhysicalVolume* GetWorldVolume();

    /**
     * Sets the safety value for the given position.
     *  @param[in] val The safety value.
     *  @param[in] pos The position.
     */
    inline void SetCurrentSafety(G4double val, const G4ThreeVector& pos);

    /**
     * Initialises all data and navigator.
     */
    void InitialiseHelper();

  private:

    G4PathFinder* fpPathFinder = nullptr;
    G4Navigator* fpMassNavigator = nullptr;

    /** Flag whether to use PathFinder or single (mass) navigator directly.
        By default, one geometry only. */
    G4bool fUseParallelGeometries = false; 

    /** Flag of first call. */
    G4bool fFirstCall = true;

    /** Whether to print warning in case of move outside safety. */
    G4int fVerbose = 0; 

    // State used during tracking -- for optimisation -------------------------

    G4ThreeVector fLastSafetyPosition;
    G4double fLastSafety = 0.0;

    // End State (tracking) ---------------------------------------------------
};

// --------------------------------------------------------------------
// Inline definitions
// --------------------------------------------------------------------

#include "G4SafetyHelper.icc"

#endif
