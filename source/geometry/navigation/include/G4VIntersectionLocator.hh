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
// G4VIntersectionLocator 
//
// Class description:
//
// Base class for the calculation of the intersection point with a boundary 
// when PropagationInField is used.
// Gives possibility to choose the method of intersection; concrete locators
// implemented are: G4SimpleLocator, G4MultiLevelLocator, G4BrentLocator.
// Key Method: EstimateIntersectionPoint().

// Authors: John Apostolakis, Tatiana Nikitina (CERN), 27 October 2008
// --------------------------------------------------------------------
#ifndef G4VINTERSECTIONLOCATOR_HH
#define G4VINTERSECTIONLOCATOR_HH 1

#include "G4Types.hh" 
#include "G4ThreeVector.hh"
#include "G4FieldTrack.hh"

#include "G4Navigator.hh"
#include "G4ChordFinder.hh"

/**
 * @brief G4VIntersectionLocator is a base class for the calculation of the
 * intersection point with a boundary when PropagationInField is used.
 * It gives the possibility to choose the method of intersection.
 */

class G4VIntersectionLocator
{
  public:
 
    /**
     * Constructor and virtual Destructor.
     */
    G4VIntersectionLocator(G4Navigator *theNavigator);
    virtual ~G4VIntersectionLocator();
     
    /**
     * If such an intersection exists, this method calculates the intersection
     * point of the true path of the particle with the surface of the current
     * volume (or of one of its daughters). 
     * Should use lateral displacement as measure of convergence.
     *  @note Changes the safety!
     *  @param[in] curveStartPointTangent Start point tangent track.
     *  @param[in] curveEndPointTangent End point tangent track.
     *  @param[in] trialPoint Trial point.
     *  @param[out] intersectPointTangent Intersection point tangent track.
     *  @param[out] recalculatedEndPoint Flagging if end point was recomputed.
     *  @param[in,out] fPreviousSafety Previous safety distance.
     *  @param[in,out] fPreviousSftOrigin Previous safety point origin.
     *  @returns Whether intersection exists or not. 
     */
    virtual G4bool EstimateIntersectionPoint( 
                     const G4FieldTrack&  curveStartPointTangent,  // A
                     const G4FieldTrack&  curveEndPointTangent,    // B
                     const G4ThreeVector& trialPoint,              // E
                           G4FieldTrack&  intersectPointTangent,   // Output
                           G4bool&        recalculatedEndPoint,    // Out
                           G4double&      fPreviousSafety,         // In/Out
                           G4ThreeVector& fPreviousSftOrigin) = 0; // In/Out   

    /**
     * Intersects the chord from StartPointA to EndPointB and returns
     * whether an intersection occurred.
     *  @note Changes the safety!
     *  @param[in] StartPointA Chord starting point.
     *  @param[in] EndPointB Chord end point.
     *  @param[out] NewSafety New calculated safety distance.
     *  @param[in,out] PreviousSafety Previous safety distance.
     *  @param[in,out] PreviousSftOrigin Previous safety point origin.
     *  @param[out] LinearStepLength Linear chord length.
     *  @param[out] IntersectionPoint Intersection point.
     *  @param[in,out] calledNavigator Pointer to flag indicating if the
     *                 navigator has been called or not.
     *  @returns Whether intersection exists or not. 
     */
    inline G4bool IntersectChord( const G4ThreeVector&  StartPointA,
                                  const G4ThreeVector&  EndPointB,
                                  G4double&      NewSafety,
                                  G4double&      PreviousSafety,    // In/Out
                                  G4ThreeVector& PreviousSftOrigin, // In/Out
                                  G4double&      LinearStepLength,
                                  G4ThreeVector& IntersectionPoint,
                                  G4bool*        calledNavigator = nullptr );

    /**
     * Setters for parameters which must be set at each step, in case they are
     * changed.
     *  @note This simple approach ensures that all scenarios are considered. 
     *        Future refinement may identify which are invariant during a 
     *        track, run or event.
     */
    inline void SetEpsilonStepFor( G4double EpsilonStep );
    inline void SetDeltaIntersectionFor( G4double deltaIntersection );
    inline void SetNavigatorFor( G4Navigator* fNavigator );
    inline void SetChordFinderFor(G4ChordFinder* fCFinder );

    /**
     * Verbosity control.
     * Controlling verbosity enables checking of the locating of intersections.
     */
    inline void  SetVerboseFor(G4int fVerbose);
    inline G4int GetVerboseFor();

    /**
     * Additional inline Get/Set methods for parameters, dependent objects.
     */
    inline G4double       GetDeltaIntersectionFor();
    inline G4double       GetEpsilonStepFor();
    inline G4Navigator*   GetNavigatorFor();
    inline G4ChordFinder* GetChordFinderFor();
    inline void SetSafetyParametersFor(G4bool UseSafety );

    /**
     * Adjustment flag accessor/modifier.
     */
    inline void AdjustIntersections(G4bool UseCorrection); 
    inline G4bool AreIntersectionsAdjusted(){ return fUseNormalCorrection; }  

    /**
     * Adjustment flag accessor/modifier.
     *  @deprecated Replaced by methods above. To be removed in future releases.
     */
    inline void AddAdjustementOfFoundIntersection(G4bool UseCorrection);
    inline G4bool GetAdjustementOfFoundIntersection();

    /**
     * Dumps status of propagator to any ostream.
     */
    static void printStatus( const G4FieldTrack& startFT,
                             const G4FieldTrack& currentFT, 
                                   G4double      requestStep, 
                                   G4double      safety,
                                   G4int         stepNum,
                                   std::ostream& oss,
                                   G4int         verboseLevel );

    /**
     * Dumps status of propagator to cout, useful mostly for debugging.
     */
    void printStatus( const G4FieldTrack& startFT,
                      const G4FieldTrack& currentFT, 
                            G4double      requestStep, 
                            G4double      safety,
                            G4int         stepNum);

    /**
     * Sets/gets check mode.
     * When enabled, uses additional verifications and stricter condictions
     * for ensuring correctness. Effective only when G4VERBOSE is enabled.
     */
    inline void SetCheckMode( G4bool value ) { fCheckMode = value; }
    inline G4bool GetCheckMode()             { return fCheckMode; }

  protected:

    /**
     * Returns new estimate for state after curveDist starting from
     * CurrentStateA, to replace EstimtdEndStateB, and reports displacement
     * (if field is compiled verbose).
     *  @param[in] CurrentStateA Start point tangent track.
     *  @param[in] EstimtdEndStateB Estimated end point tangent track.
     *  @param[in] linearDistSq Not used.
     *  @param[in] curveDist Not used.
     *  @returns New estimate for state.
     */
    G4FieldTrack ReEstimateEndpoint( const G4FieldTrack& CurrentStateA,  
                                     const G4FieldTrack& EstimtdEndStateB,
                                           G4double linearDistSq, // not used
                                           G4double curveDist );  // not used 

    /**
     * Checks whether EndB is too far from StartA to be reached and if,
     * re-estimates new value for EndB (return in RevisedEndPoint).
     * Reports error if EndB is before StartA (in curve length)
     * In that case return errorCode = 2.
     *  @param[in] CurrentStartA Start point tangent track.
     *  @param[in] EstimatedEndB Estimated end point tangent track.
     *  @param[in] RevisedEndPoint Revised end point tangent track.
     *  @param[in] errorCode Error code (0=OK, 1=coincident points, 2=error).
     *  @returns true if end point has been revised.
     */
    G4bool CheckAndReEstimateEndpoint( const G4FieldTrack& CurrentStartA,  
                                       const G4FieldTrack& EstimatedEndB,
                                             G4FieldTrack& RevisedEndPoint,
                                             G4int&        errorCode);

    /**
     * Returns the surface normal. Position *must* be the intersection point
     * from last call to G4Navigator's ComputeStep() (via IntersectChord).
     * It tries to use cached (last) value in Navigator for speed, if it was
     * kept and valid. The value returned is in global coordinates.
     *  @note It does NOT guarantee to obtain Normal. This can happen e.g.
     *        if the "Intersection" Point is not on a surface, potentially
     *        due to either inaccuracies in the transformations used, or
     *        issues with the Solid.
     *  @param[in] CurrentInt_Point Current point.
     *  @param[in,out] validNormal Flagging if normal is a valid vector.
     *  @returns The surface normal vector in local coordinates.
     */
    G4ThreeVector GetSurfaceNormal(const G4ThreeVector& CurrentInt_Point,
                                         G4bool& validNormal);

    /**
     * Returns the surface normal of the Intersecting Solid in global
     * coordinates.
     *  @note This method is costlier then GetSurfaceNormal().
     *  @param[in] CurrentInt_Point Current point.
     *  @param[in,out] validNormal Flagging if normal is a valid vector.
     *  @returns The surface normal vector in global coordinates.
     */
    G4ThreeVector GetGlobalSurfaceNormal(const G4ThreeVector& CurrentE_Point,
                                               G4bool& validNormal);
    /**
     * Optional method for adjustment of located intersection point using
     * the surface-normal.
     *  @param[in] A Chord starting point.
     *  @param[in] CurrentE_Point E Chord point.
     *  @param[in] CurrentF_Point F Chord point.
     *  @param[in] MomentumDir Momentum direction.
     *  @param[in] IntersectAF First part intersecting?
     *  @param[in,out] IntersectionPoint Intersection point tangent track.
     *  @param[in,out] NewSafety New safety distance.
     *  @param[in,out] fPrevSafety Previous safety distance.
     *  @param[in,out] fPrevSftOrigin Previous safety point origin.
     *  @returns Whether intersection exists or not. 
     */
    G4bool AdjustmentOfFoundIntersection(const G4ThreeVector& A,
                                         const G4ThreeVector& CurrentE_Point, 
                                         const G4ThreeVector& CurrentF_Point,
                                         const G4ThreeVector& MomentumDir,
                                         const G4bool         IntersectAF, 
                                               G4ThreeVector& IntersectionPoint,
                                               G4double&      NewSafety,
                                               G4double&      fPrevSafety,
                                               G4ThreeVector& fPrevSftOrigin );
  
    /**
     * Prints a three-line report on the current "sub-step",
     * i.e. trial intersection.
     *  @param[in] step_no Step number.
     *  @param[in] ChordAB_v AB chord.
     *  @param[in] ChordEF_v EF chord.
     *  @param[in] NewMomentumDir Momentum direction.
     *  @param[in] NormalAtEntry Normal vector at entry.
     *  @param[in] validNormal Validity flag for normal vector at E.
     */
    void ReportTrialStep( G4int step_no, 
                          const G4ThreeVector& ChordAB_v,
                          const G4ThreeVector& ChordEF_v,
                          const G4ThreeVector& NewMomentumDir,
                          const G4ThreeVector& NormalAtEntry,
                          G4bool validNormal   );

    /**
     * Locates a point using the navigator and updates the state of Navigator.
     * By default, it assumes that the point is inside the current volume,
     * and returns true.
     * In check mode, it checks whether the point is *inside* the volume.
     * If it is inside, it returns true.
     * If not, issues a warning and returns false.
     *  @param[in] pos The point to locate.
     *  @returns If a point is inside the volume or not.
     */
    G4bool LocateGlobalPointWithinVolumeAndCheck( const G4ThreeVector& pos );

    /**
     * Locates a point using the navigator and updates the state of Navigator,
     * but report information about code location.
     * If CheckMode > 1, report extra information.
     *  @param[in] pos The point to locate.
     *  @param[in] CodeLocationInfo String for code location info.
     *  @param[in] CheckMode Not used.
     */
    void LocateGlobalPointWithinVolumeCheckAndReport( const G4ThreeVector& pos,
                                            const G4String& CodeLocationInfo,
                                                  G4int     CheckMode );

    // Auxiliary methods -- to report issues

    /**
     * Builds error message (in ossMsg) to report that point 'B' has
     * gone past 'A'.
     */
    void ReportReversedPoints( std::ostringstream& ossMsg,
                               const G4FieldTrack& StartPointVel, 
                               const G4FieldTrack& EndPointVel,
                                     G4double NewSafety, G4double epsStep,
                               const G4FieldTrack& CurrentA_PointVelocity,
                               const G4FieldTrack& CurrentB_PointVelocity,
                               const G4FieldTrack& SubStart_PointVelocity,
                               const G4ThreeVector& CurrentE_Point,
                               const G4FieldTrack& ApproxIntersecPointV,
                               G4int sbstp_no, G4int sbstp_no_p, G4int depth );

    /**
     * Reports the current status / progress in finding the first intersection.
     */
    void ReportProgress( std::ostream& oss,
                         const G4FieldTrack& StartPointVel, 
                         const G4FieldTrack& EndPointVel,
                               G4int         substep_no, 
                         const G4FieldTrack& A_PtVel,    // G4double safetyA
                         const G4FieldTrack& B_PtVel,  
                               G4double      safetyLast,
                               G4int         depth= -1 );

    /**
     * Report case: trial point is 'close' to start, within tolerance.
     */
     void ReportImmediateHit( const char*          MethodName, 
                              const G4ThreeVector& StartPosition, 
                              const G4ThreeVector& TrialPoint, 
                                    G4double       tolerance,
                                 unsigned long int numCalls );
    
  private:

    /**
     * Returns the SurfaceNormal of the Intersecting Solid in local coordinates.
     */
    G4ThreeVector GetLocalSurfaceNormal(const G4ThreeVector& CurrentE_Point,
                                              G4bool& validNormal);

    /**
     * Position *must* be the intersection point from last call
     * to G4Navigator's ComputeStep  (via IntersectChord).
     */
    G4ThreeVector GetLastSurfaceNormal( const G4ThreeVector& intersectPoint,
                                              G4bool& validNormal) const;

  protected:

    G4double kCarTolerance;                  // Constant

    G4int    fVerboseLevel = 0;              // For debugging
    G4bool   fUseNormalCorrection = false;   // Configuration parameter
    G4bool   fCheckMode = false;
    G4bool   fiUseSafety = false;    // Whether to use safety for 'fast steps'
   
    G4Navigator* fiNavigator;

    /**
     * Parameters set at each physical step by G4PropagatorInField.
     */
    G4ChordFinder* fiChordFinder = nullptr;  // Overridden at each step
    G4double fiEpsilonStep = -1.0;           // Overridden at each step
    G4double fiDeltaIntersection = -1.0;     // Overridden at each step

    /** Helper for location. */
    G4Navigator *fHelpingNavigator;

    /** Touchable history hook. */
    G4TouchableHistory *fpTouchable = nullptr;
};

#include "G4VIntersectionLocator.icc"

#endif
