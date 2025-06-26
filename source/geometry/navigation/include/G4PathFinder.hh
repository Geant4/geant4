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
// G4PathFinder 
//
// Class description:
// 
// This class directs the lock-stepped propagation of a track in the 
// 'mass' and other parallel geometries. It ensures that tracking 
// in a magnetic field sees these parallel geometries at each trial step, 
// and that the earliest boundary limits the step.
// 
// In field, it relies on the class G4PropagatorInField.

// Author: John Apostolakis (CERN), 7 October 2005
// ---------------------------------------------------------------------
#ifndef G4PATHFINDER_HH 
#define G4PATHFINDER_HH  1

#include <vector>

#include "G4Types.hh"
#include "G4FieldTrack.hh"
#include "G4MultiNavigator.hh"
#include "G4TouchableHandle.hh"

class G4TransportationManager; 
class G4Navigator;
class G4PropagatorInField;

/**
 * @brief G4PathFinder directs the lock-stepped propagation of a track in the 
 * 'mass' and other parallel geometries. It ensures that tracking in a magnetic
 * field sees these parallel geometries at each trial step and that the earliest
 * boundary limits the step.
 */

class G4PathFinder
{
  public:

    /**
     * Retrieves singleton instance and creates it if not existing.
     */
    static G4PathFinder* GetInstance();

    /**
     * Retrieve singleton instance pointer.
     */
    static G4PathFinder* GetInstanceIfExist();

    /**
     * Destructor, called only by G4RunManagerKernel.
     */
    ~G4PathFinder();

    /**
     * Computes the next geometric Step, curved or linear.
     * If it is called with a larger 'stepNo' it will execute a new step;
     * if 'stepNo' is same as last call, then the results for the geometry
     * with Id number 'navigatorId' will be returned.
     *  @param[in,out] pFieldTrack Field track to be filled.
     *  @param[in] pCurrentProposedStepLength Current proposed step length.
     *  @param[in] navigatorId Identifier of the geometry.
     *  @param[in] stepNo Step number; see next step/check.
     *  @param[in, out] pNewSafety New safety for this geometry.
     *  @param[in, out] limitedStep Step characterisation to be returned.
     *  @param[in, out] EndState Field track end state.
     *  @param[in] currentVolume Pointer to the current volume.
     *  @returns Step length.
     */
    G4double ComputeStep( const G4FieldTrack& pFieldTrack,
                          G4double  pCurrentProposedStepLength,
                          G4int     navigatorId, // Identifies the geometry
                          G4int     stepNo,      // See next step/check
                          G4double& pNewSafety,  // Only for this geometry
                          ELimited& limitedStep,
                          G4FieldTrack& EndState, 
                          G4VPhysicalVolume* currentVolume );

    /**
     * Makes primary relocation of the global point in all navigators,
     * and updates them.
     *  @param[in] position Point in global coordinates system.
     *  @param[in] direction Global direction vector.
     *  @param[in] relativeSearch If set to true (default), the search begins
     *             is the geometrical hierarchy at the location of the last
     *             located point.
     */
    void Locate( const G4ThreeVector& position, 
                 const G4ThreeVector& direction,
                       G4bool relativeSearch = true); 

    /**
     * Makes secondary relocation of the global point (within safety only) 
     * in all navigators, and updates them.
     *  @param[in] position Point in global coordinates system.
     */
    void ReLocate( const G4ThreeVector& position ); 

    /**
     * Checks and caches the set of active navigators.
     *  @param[in] position Point in global coordinates system.
     *  @param[in] direction Global direction vector.
     *  @param[in] massStartVol Pointer to the mass geometry world.
     */
    void PrepareNewTrack( const G4ThreeVector& position,
                          const G4ThreeVector& direction,
                                G4VPhysicalVolume* massStartVol = nullptr);

    /**
     * Signals the end of tracking of the current track. Resets internal state
     * and informs G4TransportationManager to use 'ordinary' Navigator.
     */
    void EndTrack();

    /**
     * Creates a touchable handle for the specified navigator.
     *  @param[in] navId The navigator identifier.
     *  @returns A touchable handle of the geometry.
     */
    G4TouchableHandle CreateTouchableHandle( G4int navId ) const;

    /**
     * Returns the located volume for the specified navigator.
     *  @param[in] navId The navigator identifier.
     *  @returns A pointer to the located volume in the geometry.
     */
    inline G4VPhysicalVolume* GetLocatedVolume( G4int navId ) const;

    // -----------------------------------------------------------------
  
    /**
     * Returns the minimum value of safety after last ComputeStep().
     */
    inline G4double GetCurrentSafety() const;

    /**
     * Gets the minimum step size from the last ComputeStep() call.
     *  @note In case full step is taken, this is kInfinity.
     */
    inline G4double GetMinimumStep() const;      

    /**
     * Returns the number of all geometries limiting the step.
     */
    inline unsigned int GetNumberGeometriesLimitingStep() const; 

    /**
     * Recomputes the safety for the relevant point, i.e. the endpoint of the
     * last step. Maintains a vector of individual safety values (used by next
     * method below).
     *  @param[in] globalPoint Point in global coordinates system.
     *  @returns The safety value for the specified point in the geometry.
     */
    G4double ComputeSafety( const G4ThreeVector& globalPoint); 

    /**
     * Obtains the safety for the specified navigator/geometry for last point
     * 'computed' (i.e., the last point for which ComputeSafety() was called).
     *  @param[in] navId The navigator identifier.
     *  @param[in,out] globalCenterPoint The point (center) for which this
     *                 safety is valid.
     *  @returns The safety value in the specified geometry.
     */
    inline G4double ObtainSafety(G4int navId, G4ThreeVector& globalCenterPoint);

    /**
     * To enable parallel navigation. Must call it to ensure that G4PathFinder
     * is prepared, especially for curved tracks.
     * If true it switches G4PropagatorInField to use G4MultiNavigator.
     * Must call it with false to undo (i.e. G4PropagatorInField uses classic
     * G4Navigator for tracking in such case).
     *  @param[in] enableChoice Flag to enable/disable parallel navigation.
     */
    void EnableParallelNavigation( G4bool enableChoice = true ); 

    /**
     * To control the level of verbosity. Default is no verbosity.
     */
    inline G4int SetVerboseLevel(G4int lev = -1);

    /**
     * To get/set the maximum for the number of steps that a (looping)
     * particle can take.
     */
    inline G4int GetMaxLoopCount() const;
    inline void  SetMaxLoopCount( G4int new_max );

    /**
     * Signals that the point location will be moved.
     *  @note Internal use primarily.
     */
    inline void MovePoint();

    // To provide best compatibility between Coupled and normal Transportation
    // the next two methods are provided...

    /**
     * Obtains the last safety needed in ComputeStep() for the specified
     * geometry 'navId' (i.e. the last point at which ComputeStep() has
     * recalculated the safety). Returns the point (center) for which this
     * safety is valid and also the minimum safety over all navigators.
     *  @param[in] navId The navigator identifier.
     *  @param[in,out] globCenterPoint The point (center) for which this
     *                 safety is valid.
     *  @param[in,out] minSafety The minimum safety over all navigators.
     *  @returns The safety value in the specified geometry.
     */
    inline G4double LastPreSafety( G4int navId, G4ThreeVector& globCenterPoint,
                                   G4double& minSafety ); 

    /**
     * Tells G4PathFinder to copy PostStep Safety to PreSafety
     * for use at the next step.
     */
    void PushPostSafetyToPreSafety(); 

    /**
     * Utility to convert ELimited specification to a string.
     */
    G4String& LimitedString( ELimited lim );

  private:

    /**
     * Private singleton constructor.
     */
    G4PathFinder();

    /**
     * Returns pointer to the specified navigator.
     */
    inline G4Navigator* GetNavigator(G4int n) const;

    /**
     * Performs a linear step.
     *  @param[in,out] FieldTrack Field track to be filled.
     *  @param[in] proposedStepLength Current proposed step length.
     *  @returns The minimum linear step to undertake.
     */
    G4double DoNextLinearStep( const G4FieldTrack& FieldTrack,
                                     G4double proposedStepLength); 

    /**
     * Performs a curved step.
     *  @param[in,out] FieldTrack Field track to be filled.
     *  @param[in] proposedStepLength Current proposed step length.
     *  @param[in] pCurrentPhysVolume Pointer to the current volume of interest.
     *  @returns The minimum step to undertake.
     */
    G4double DoNextCurvedStep( const G4FieldTrack& FieldTrack,
                                     G4double proposedStepLength,
                                     G4VPhysicalVolume* pCurrentPhysVolume); 

    /**
     * Prints key details out for debugging.
     */
    void WhichLimited();
    void PrintLimited();

    /**
     * Helper method to report movement (likely of initial point).
     */
    void ReportMove( const G4ThreeVector& OldV,
                     const G4ThreeVector& NewV,
                     const G4String& Quantity ) const;

  private:

    // ----------------------------------------------------------------------
    //  DATA Members
    // ----------------------------------------------------------------------

    /** Object that enables G4PropagatorInField to see many geometries. */
    G4MultiNavigator* fpMultiNavigator; 

    G4int fNoActiveNavigators = 0; 
    G4bool fNewTrack = false; // Flag a new track (ensure first step)

    static const G4int fMaxNav = 16;

    // Global state (retained during stepping for one track)

    G4Navigator*  fpNavigator[fMaxNav];

    // ---- State changed in a step computation
    //
    ELimited fLimitedStep[fMaxNav];
    G4bool fLimitTruth[fMaxNav];
    G4double fCurrentStepSize[fMaxNav]; 
    G4int fNoGeometriesLimiting = 0;  // How many processes contribute to limit

    /** Last initial position for which safety evaluated. */
    G4ThreeVector fPreSafetyLocation;
    /* Corresponding value of full safety. */
    G4double fPreSafetyMinValue = -1.0;

    /** Safeties for the above point. */
    G4double fPreSafetyValues[ fMaxNav ];

    // This part of the state can be retained for several calls --> CARE

    /** Point where last ComputeStep() called. */
    G4ThreeVector fPreStepLocation;
    /** Corresponding value of full safety. */
    G4double fMinSafety_PreStepPt = -1.0;

    /** Safeties for the above point.
      * @note This changes at each step, so it can differ when steps
      *       inside min-safety are made. */
    G4double fCurrentPreStepSafety[ fMaxNav ];

    /** Whether PreSafety coincides with PreStep point. */
    G4bool fPreStepCenterRenewed = false;

    G4double fMinStep = -1.0;  // As reported by Navigators -- can be kInfinity
    G4double fTrueMinStep = -1.0;  // Corrected in case >= proposed

    // ---- State after calling 'locate'
    //
    G4VPhysicalVolume* fLocatedVolume[fMaxNav];
    G4ThreeVector      fLastLocatedPosition; 

    // ---- State after calling 'ComputeStep'
    //      (other member variables will be affected)
    //
    G4FieldTrack fEndState;  // Point, velocity, ... at proposed step end
    G4bool fFieldExertedForce = false; // In current proposed step

    G4bool fRelocatedPoint = false; // Signals that point was or is being moved 
                                    // from the position of the last location or
                                    // the endpoint resulting from ComputeStep()
                                    // -- invalidates fEndState

    // ---- State for 'ComputeSafety' and related methods
    //
    /** Point where ComputeSafety() is called. */
    G4ThreeVector fSafetyLocation;
    /** Corresponding value of safety. */
    G4double      fMinSafety_atSafLocation = -1.0;
    /** Safeties for last ComputeSafety(). */
    G4double      fNewSafetyComputed[ fMaxNav ];

    // ---- State for Step numbers
    //
    G4int fLastStepNo = -1, fCurrentStepNo = -1; 

    G4int fVerboseLevel = 0;  // For debugging purposes

    G4TransportationManager* fpTransportManager; // Cache for frequent use
    G4PropagatorInField* fpFieldPropagator;

    G4double kCarTolerance;

    static G4ThreadLocal G4PathFinder* fpPathFinder;
};

// ********************************************************************
// Inline methods.
// ********************************************************************

#include "G4PathFinder.icc"

#endif 
