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
// $Id: G4VIntersectionLocator.hh 93289 2015-10-15 10:01:15Z gcosmo $
//
//
// Class G4VIntersectionLocator 
//
// class description:
//
// Base class for the calculation of the intersection point with a boundary 
// when PropagationInField is used.
// Gives possibility to choose the method of intersection; concrete locators
// implemented are: G4SimpleLocator, G4MultiLevelLocator, G4BrentLocator.
//
// Key Method: EstimateIntersectionPoint()

// History:
// -------
// 27.10.08 - John Apostolakis, Tatiana Nikitina: Design and implementation 
// ---------------------------------------------------------------------------

#ifndef G4VINTERSECTIONLOCATOR_HH
#define G4VINTERSECTIONLOCATOR_HH

#include "G4Types.hh" 
#include "G4ThreeVector.hh"
#include "G4FieldTrack.hh"

#include "G4Navigator.hh"
#include "G4ChordFinder.hh"

class G4VIntersectionLocator
 {
   public:  // with description 
 
     G4VIntersectionLocator(G4Navigator *theNavigator);
       // Constructor
     virtual ~G4VIntersectionLocator();
       // Default destructor
     
     virtual G4bool EstimateIntersectionPoint( 
         const  G4FieldTrack&       curveStartPointTangent,  // A
         const  G4FieldTrack&       curveEndPointTangent,    // B
         const  G4ThreeVector&      trialPoint,              // E
                G4FieldTrack&       intersectPointTangent,   // Output
                G4bool&             recalculatedEndPoint,    // Out
                G4double&           fPreviousSafety,         // In/Out
                G4ThreeVector&      fPreviousSftOrigin) = 0; // In/Out   
       // If such an intersection exists, this function calculates the
       // intersection point of the true path of the particle with the surface
       // of the current volume (or of one of its daughters). 
       // Should use lateral displacement as measure of convergence
       // NOTE: changes the safety!

     void printStatus( const G4FieldTrack& startFT,
                       const G4FieldTrack& currentFT, 
                             G4double      requestStep, 
                             G4double      safety,
                             G4int         stepNum);
       // Print Method, useful mostly for debugging

     inline G4bool IntersectChord( const G4ThreeVector&  StartPointA,
                                   const G4ThreeVector&  EndPointB,
                                   G4double      &NewSafety,
                                   G4double      &PreviousSafety,    // In/Out
                                   G4ThreeVector &PreviousSftOrigin, // In/Out
                                   G4double      &LinearStepLength,
                                   G4ThreeVector &IntersectionPoint,
                                   G4bool        *calledNavigator=0 );
       // Intersect the chord from StartPointA to EndPointB and return
       // whether an intersection occurred. NOTE: changes the Safety!

     inline void    SetEpsilonStepFor( G4double EpsilonStep );
     inline void    SetDeltaIntersectionFor( G4double deltaIntersection );
     inline void    SetNavigatorFor( G4Navigator *fNavigator );
     inline void    SetChordFinderFor(G4ChordFinder *fCFinder );
       // These parameters must be set at each step, in case they were changed

       // Note: This simple approach ensures that all scenarios are considered. 
       //   [ Future refinement may identify which are invariant during a 
       //      track, run or event ]

    inline void     SetVerboseFor(G4int fVerbose);
    inline G4int    GetVerboseFor();
       // Controling verbosity enables checking of the locating of intersections

  public:  // without description

    // Additional inline Set/Get methods for parameters, dependent objects

    inline G4double       GetDeltaIntersectionFor();
    inline G4double       GetEpsilonStepFor();
    inline G4Navigator*   GetNavigatorFor();
    inline G4ChordFinder* GetChordFinderFor();

    inline void   SetSafetyParametersFor(G4bool UseSafety );

    inline void   AddAdjustementOfFoundIntersection(G4bool UseCorrection);
    inline G4bool GetAdjustementOfFoundIntersection();
      // Methods to be made Obsolete - replaced by methods below
    inline void   AdjustIntersections(G4bool UseCorrection); 
    inline G4bool AreIntersectionsAdjusted(){ return fUseNormalCorrection; }  
      // Change adjustment flag  ( New Interface ) 

    static void printStatus( const G4FieldTrack& startFT,
                             const G4FieldTrack& currentFT, 
                                   G4double      requestStep, 
                                   G4double      safety,
                                   G4int         stepNum,
                                   std::ostream& oss,
                                   G4int         verboseLevel );
      // Print Method for any ostream - e.g. cerr -- and for G4Exception
    
  protected:  // with description

    G4FieldTrack ReEstimateEndpoint( const G4FieldTrack &CurrentStateA,  
                                     const G4FieldTrack &EstimtdEndStateB,
                                           G4double linearDistSq, // not used
                                           G4double curveDist );  // not used
      // Return new estimate for state after curveDist starting from
      // CurrentStateA, to replace EstimtdEndStateB, and report displacement
      // (if field is compiled verbose)

    G4bool CheckAndReEstimateEndpoint( const G4FieldTrack& CurrentStartA,  
                                       const G4FieldTrack& EstimatedEndB,
                                             G4FieldTrack& RevisedEndPoint,
                                             G4int &       errorCode);
      // Check whether EndB is too far from StartA to be reached 
      // and if, re-estimate new value for EndB (return in RevisedEndPoint)
      // Report error if  EndB is before StartA (in curve length)
      // In that case return errorCode = 2.

    G4ThreeVector GetSurfaceNormal(const G4ThreeVector &CurrentInt_Point,
                                         G4bool        &validNormal); // const
      // Position *must* be the intersection point from last call
      // to G4Navigator's ComputeStep  (via IntersectChord )
      // Will try to use cached (last) value in Navigator for speed, 
      // if it was kept and valid.
      // Value returned is in global coordinates.
      // It does NOT guarantee to obtain Normal. This can happen eg if:
      //  - the "Intersection" Point is not on a surface, potentially due to
      //  - inaccuracies in the transformations used, or
      //  - issues with the Solid.

    G4ThreeVector GetGlobalSurfaceNormal(const G4ThreeVector &CurrentE_Point,
                                               G4bool        &validNormal);
      // Return the SurfaceNormal of Intersecting Solid in global coordinates
      // Costlier then GetSurfaceNormal

    G4bool AdjustmentOfFoundIntersection(const G4ThreeVector &A,
                                         const G4ThreeVector &CurrentE_Point, 
                                         const G4ThreeVector &CurrentF_Point,
                                         const G4ThreeVector &MomentumDir,
                                         const G4bool         IntersectAF, 
                                               G4ThreeVector &IntersectionPoint,
                                               G4double      &NewSafety,
                                               G4double      &fPrevSafety,
                                               G4ThreeVector &fPrevSftOrigin );
      // Optional method for adjustment of located intersection point
      // using the surface-normal
  
    void ReportTrialStep( G4int step_no, 
                          const G4ThreeVector& ChordAB_v,
                          const G4ThreeVector& ChordEF_v,
                          const G4ThreeVector& NewMomentumDir,
                          const G4ThreeVector& NormalAtEntry,
                          G4bool validNormal   );
      // Print a three-line report on the current "sub-step", ie trial
      // intersection

    G4bool LocateGlobalPointWithinVolumeAndCheck( const G4ThreeVector& pos );
      // Locate point using navigator - updates state of Navigator.
      // By default, it assumes that the point is inside the current volume,
      // and returns true.
      // In check mode, checks whether the point is *inside* the volume.
      //   If it is inside, it returns true.
      //   If not, issues a warning and returns false.

    void LocateGlobalPointWithinVolumeCheckAndReport( const G4ThreeVector& pos,
                                            const G4String& CodeLocationInfo,
                                                  G4int     CheckMode );
      // As above, but report information about code location.
      // If CheckMode > 1, report extra information.
   
    inline void   SetCheckMode( G4bool value ) { fCheckMode = value; }
    inline G4bool GetCheckMode()               { return fCheckMode; }

  protected:  // without description

    // Auxiliary methods -- to report issues

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
      // Build error messsage (in ossMsg) to report that point 'B' has
      // gone past 'A'

    void ReportProgress( std::ostream& oss,
                         const G4FieldTrack& StartPointVel, 
                         const G4FieldTrack& EndPointVel,
                               G4int         substep_no, 
                         const G4FieldTrack& A_PtVel,    // G4double safetyA
                         const G4FieldTrack& B_PtVel,  
                               G4double      safetyLast,
                               G4int         depth= -1 );
      // Report the current status / progress in finding the first intersection

     void ReportImmediateHit( const char*          MethodName, 
                              const G4ThreeVector& StartPosition, 
                              const G4ThreeVector& TrialPoint, 
                                    double         tolerance,
                                 unsigned long int numCalls );
      // Report case: trial point is 'close' to start, within tolerance
    
  private:  // no description

    G4ThreeVector GetLocalSurfaceNormal(const G4ThreeVector &CurrentE_Point,
                                              G4bool &validNormal);
      // Return the SurfaceNormal of Intersecting Solid  in local coordinates

    G4ThreeVector GetLastSurfaceNormal( const G4ThreeVector& intersectPoint,
                                        G4bool        &validNormal) const; 
      // Position *must* be the intersection point from last call
      // to G4Navigator's ComputeStep  (via IntersectChord )
      // Temporary - will use the same method in the Navigator

  protected:

    G4double kCarTolerance;         // Constant

    G4int    fVerboseLevel;          // For debugging
    G4bool   fUseNormalCorrection;   // Configuration parameter
    G4bool   fCheckMode;
   
    G4Navigator   *fiNavigator;

    G4ChordFinder *fiChordFinder;
    G4double       fiEpsilonStep;
    G4double       fiDeltaIntersection;
    G4bool         fiUseSafety;
      // Parameters set at each physical step by calling method 
      // by G4PropagatorInField

    G4Navigator *fHelpingNavigator;
      // Helper for location

    G4TouchableHistory *fpTouchable;
      // Touchable history hook
};

#include "G4VIntersectionLocator.icc"

#endif
