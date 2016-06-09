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
// $Id: G4PropagatorInField.hh,v 1.13 2007/06/08 09:49:34 gcosmo Exp $
// GEANT4 tag $Name: geant4-09-00 $
// 
// class G4PropagatorInField 
//
// Class description:
// 
// This class performs the navigation/propagation of a particle/track 
// in a magnetic field. The field is in general non-uniform.
// For the calculation of the path, it relies on the class G4ChordFinder.
//
// Key Method:
//              ComputeStep(..)
// History:
// -------
// 25.10.96 John Apostolakis,  design and implementation 
// 25.03.97 John Apostolakis,  adaptation for G4Transportation and cleanup
//  8.11.02 John Apostolakis,  changes to enable use of safety in intersecting
// ---------------------------------------------------------------------------

#ifndef G4PropagatorInField_hh 
#define G4PropagatorInField_hh  1

#include "G4Types.hh"

#include <vector>

#include "G4FieldTrack.hh"
#include "G4FieldManager.hh"
class G4ChordFinder; 

class G4Navigator;
class G4VPhysicalVolume;
class G4VCurvedTrajectoryFilter;

class G4PropagatorInField
{

 public:  // with description

   G4PropagatorInField( G4Navigator    *theNavigator, 
                        G4FieldManager *detectorFieldMgr );
  ~G4PropagatorInField();

   G4double ComputeStep( G4FieldTrack      &pFieldTrack,
                         G4double           pCurrentProposedStepLength,
                         G4double          &pNewSafety, 
                         G4VPhysicalVolume *pPhysVol=0 );
     // Compute the next geometric Step

   inline G4ThreeVector  EndPosition() const;       
   inline G4ThreeVector  EndMomentumDir() const;
   inline G4bool         IsParticleLooping() const;
     // Return the state after the Step

   inline G4double  GetEpsilonStep() const;
     // Relative accuracy for current Step (Calc.)
   inline void      SetEpsilonStep(G4double newEps);
     // The ratio DeltaOneStep()/h_current_step

   inline void SetChargeMomentumMass( G4double charge,     // in e+ units
                                      G4double momentum,   // in Geant4 units
                                      G4double pMass );  
     // Inform this and all associated classes of q, p, m 

   G4FieldManager*  FindAndSetFieldManager(G4VPhysicalVolume* pCurrentPhysVol);
     // Set (and return) the correct field manager (global or local), 
     //    if it exists.
     // Should be called before ComputeStep is called;
     //   - currently, ComputeStep will call it, if it has not been called.
 
   inline G4ChordFinder* GetChordFinder();

          G4int  SetVerboseLevel( G4int verbose );
   inline G4int  GetVerboseLevel() const;
   inline G4int  Verbose() const;

   inline G4int   GetMaxLoopCount() const;
   inline void    SetMaxLoopCount( G4int new_max );
     // A maximum for the number of steps that a (looping) particle can take.

   void printStatus( const G4FieldTrack&        startFT,
                     const G4FieldTrack&        currentFT, 
                           G4double             requestStep, 
                           G4double             safety,
                           G4int                step, 
                           G4VPhysicalVolume*   startVolume);
     // Print Method - useful mostly for debugging.

   inline G4FieldTrack GetEndState() const;

   // The following methods are now obsolescent but *for now* will work 
   //   They are being replaced by same-name methods in G4FieldManager,
   //   allowing the specialisation in different volumes. 
   //   Their new behaviour is to change the values for the global field manager. 
   inline G4double  GetMinimumEpsilonStep() const;
   inline void      SetMinimumEpsilonStep( G4double newEpsMin );
     // Minimum for Relative accuracy of any Step 

   inline G4double  GetMaximumEpsilonStep() const;
   inline void      SetMaximumEpsilonStep( G4double newEpsMax );

   inline void      SetLargestAcceptableStep( G4double newBigDist );
   inline G4double  GetLargestAcceptableStep();

 public:  // with description

   // The following methods are obsolete and will not work --
   //   as they have been replaced by the same methods in G4FieldManager
   //   since Geant4 4.0
   inline G4double  GetDeltaIntersection() const;
   inline G4double  GetDeltaOneStep() const;
   inline void    SetAccuraciesWithDeltaOneStep( G4double deltaOneStep );  
   inline void    SetDeltaIntersection( G4double deltaIntersection );
   inline void    SetDeltaOneStep( G4double deltaOneStep );  

 public:  // without description

   inline G4FieldManager*  GetCurrentFieldManager();
   inline void             SetNavigatorForPropagating( G4Navigator *SimpleOrMultiNavigator ); 
   inline G4Navigator*     GetNavigatorForPropagating(); 

 public:  // no description

   inline void SetThresholdNoZeroStep( G4int noAct,
                                       G4int noHarsh,
                                       G4int noAbandon );
   inline G4int GetThresholdNoZeroSteps( G4int i ); 

 public:  // with description
  // 
  void SetTrajectoryFilter(G4VCurvedTrajectoryFilter* filter);
  // Set the filter that examines & stores 'intermediate' 
  //  curved trajectory points.  Currently only position is stored.

  std::vector<G4ThreeVector>* GimmeTrajectoryVectorAndForgetIt() const;
  // Access the points which have passed by the filter.
  // Responsibility for deleting the points lies with the client.
  // This method MUST BE called exactly ONCE per step. 

  void ClearPropagatorState();
  // Clear all the State of this class and its current associates
  //   --> the current field manager & chord finder will also be called

  inline void SetDetectorFieldManager( G4FieldManager* newGlobalFieldManager );
      // Update this (dangerous) state -- for the time being
  
  inline void   SetUseSafetyForOptimization( G4bool );
  inline G4bool GetUseSafetyForOptimization();
      // Toggle & view parameter for using safety to discard 
      //   unneccesary calls to navigator (thus 'optimising' performance)

 protected:  // with description

   G4bool LocateIntersectionPoint( 
        const  G4FieldTrack&       curveStartPointTangent,  //  A
        const  G4FieldTrack&       curveEndPointTangent,    //  B
        const  G4ThreeVector&      trialPoint,              //  E
               G4FieldTrack&       intersectPointTangent,   // Output
               G4bool&             recalculatedEndPoint);   // Out: 

     // If such an intersection exists, this function 
     // calculate the intersection point of the true path of the particle 
     // with the surface of the current volume (or of one of its daughters). 
     // (Should use lateral displacement as measure of convergence). 

   G4bool IntersectChord( G4ThreeVector  StartPointA, 
                          G4ThreeVector  EndPointB,
                          G4double      &NewSafety,
                          G4double      &LinearStepLength,
                          G4ThreeVector &IntersectionPoint);
     // Intersect the chord from StartPointA to EndPointB
     // and return whether an intersection occurred

   G4FieldTrack ReEstimateEndpoint( const G4FieldTrack &CurrentStateA,  
                                    const G4FieldTrack &EstimtdEndStateB,
                                          G4double      linearDistSq,
                                          G4double      curveDist);
     // Return new estimate for state after curveDist 
     // starting from CurrentStateA,  to replace EstimtdEndStateB,
     // (and report displacement -- if field is compiled verbose.)

   void PrintStepLengthDiagnostic( G4double      currentProposedStepLength,
                                   G4double      decreaseFactor,
                                   G4double      stepTrial,
                             const G4FieldTrack& aFieldTrack);
 private:

  // ----------------------------------------------------------------------
  //  DATA Members
  // ----------------------------------------------------------------------

   G4FieldManager *fDetectorFieldMgr; 
     // The  Field Manager of the whole Detector.  (default)

   G4FieldManager *fCurrentFieldMgr;
     // The  Field Manager of the current volume (may be the one above.)

   G4Navigator   *fNavigator;

   //  STATE information
   //  -----------------

   G4double    fEpsilonStep;
     // Relative accuracy for current Step (Calc.)

   G4FieldTrack    End_PointAndTangent;
     // End point storage

   G4bool      fParticleIsLooping;

   G4int  fVerboseLevel;
     // For debuging purposes

   G4int  fMax_loop_count;
     // Limit for the number of sub-steps taken in one call to ComputeStep

   //  Variables to keep track of "abnormal" case - which causes loop
   //
   G4int     fNoZeroStep;                        //  Counter of zeroStep
   G4int     fActionThreshold_NoZeroSteps;       //  Threshold: above this - act
   G4int     fSevereActionThreshold_NoZeroSteps; //  Threshold to act harshly
   G4int     fAbandonThreshold_NoZeroSteps;      //  Threshold to abandon

   G4double  fFull_CurveLen_of_LastAttempt; 
   G4double  fLast_ProposedStepLength; 
   G4double  fLargestAcceptableStep;

   G4double  fCharge, fInitialMomentumModulus, fMass;

   G4ThreeVector  fPreviousSftOrigin;
   G4double       fPreviousSafety; 
   G4bool         fUseSafetyForOptimisation;
     // Last safety origin & value: for optimisation

   G4bool fSetFieldMgr; 
     // Flag whether field manager has been set for the current step

   G4double  kCarTolerance;
     // Geometrical tolerance defining surface thickness

private:

   static const G4int max_depth=4;
   G4FieldTrack* ptrInterMedFT[max_depth+1];
     // Used to store intermediate values of tracks in case of
     // too slow progress

private:

  G4VCurvedTrajectoryFilter* fpTrajectoryFilter;
    // The filter encapsulates the algorithm which selects which
    // intermediate points should be stored in a trajectory. 
    // When it is NULL, no intermediate points will be stored.
    // Else PIF::ComputeStep must submit (all) intermediate
    // points it calculates, to this filter.  (jacek 04/11/2002)
};

// ********************************************************************
// Inline methods.
// ********************************************************************

#include "G4PropagatorInField.icc"

#endif 
