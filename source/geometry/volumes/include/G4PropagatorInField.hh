// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4PropagatorInField.hh,v 1.7 2000-04-25 16:15:03 gcosmo Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// class G4PropagatorInField 
//
// Class description:
// 
// This class performs the navigation/propagation of a particle/track 
// in a magnetic field. The field is in general non-uniform.
// For the calculation of the path, it relies on the class G4MagTr.
// To create an object, must have an object that calculates the Curved 
// paths and also must know the value of the maximum displacement allowed.
//
// Methods:
//        ComputeStep(..)
//        CalculateStepTimeAndAccuracy(..) 
//        LocateIntersectionPoint(..)

// History:
// -------
// 25.10.96 John Apostolakis,  design and implementation 
// 25.03.97 John Apostolakis,  adaptation for G4Transportation and cleanup
// ------------------------------------------------------------------------

#ifndef G4PropagatorInField_hh 
#define G4PropagatorInField_hh  1

#include "globals.hh"
#include "G4FieldTrack.hh"
// #include "G4VPhysicalVolume.hh"
// class  G4VPhysicalVolume;

#include "G4Navigator.hh"
#include "G4ChordFinder.hh"

#include "G4FieldManager.hh"

// #include "G4MagIntegratorDriver.hh"

class G4PropagatorInField
{

 public:  // with description

   G4PropagatorInField( G4Navigator    *theNavigator, 
			G4FieldManager *detectorFieldMgr);

   G4PropagatorInField( G4Navigator   *theNavigator );
  ~G4PropagatorInField(){};

   G4double ComputeStep(G4FieldTrack  &pFieldTrack,
			G4double pCurrentProposedStepLength,
			G4double       &pNewSafety, 
			G4VPhysicalVolume *pPhysVol=0 );

   G4double ComputeStep(const G4ThreeVector &pGlobalPoint,
		        const G4ThreeVector &pCurveTangent,    // Unit vector
		              G4double pCurrentProposedStepLength,
			      G4double       &pNewSafety, 
                       	      G4VPhysicalVolume *pPhysVol=0 );
     // Compute the next geometric Step

   G4ThreeVector  EndPosition();       
   G4ThreeVector  EndMomentumDir();
   G4bool         IsParticleLooping();
     // Return the state after the Step

   G4double  DeltaIntersection(); 
     // The accuracy of finding an intersection

   G4double  DeltaOneStep();
     // The accuracy of a single Step

   G4double  GetEpsilonStep();  // Relative accuracy for current Step (Calc.)
   void      SetEpsilonStep(G4double newEps);
     // The ratio DeltaOneStep()/h_current_step

   void SetChargeMomentumMass( G4double Charge,         // in e+ units
			       G4double Momentum,       // in Geant4 units
			       G4double pMass);  

   G4ChordFinder* GetChordFinder();
   // void        SetChordFinder(G4ChordFinder* newCF);  // Not yet relevant

   G4int  SetVerboseLevel( G4int Verbose );
   G4int  Verbose();

   G4double  GetDeltaIntersection();
     // Accuracy for boundary intersection.
   G4double  GetDeltaOneStep();
     // Accuracy for one tracking/physics step.
                                    
   void    SetAccuraciesWithDeltaOneStep(G4double deltaOneStep);  
     // Sets both accuracies, maintaining a particular ratio
     // Delta Interaction / OneStep.

   G4int   GetMaxLoopCount();
   void    SetMaxLoopCount(G4int new_max);
     // A maximum for the number of steps that a (looping) particle can take.

   void printStatus( 
                  const G4FieldTrack&  StartFT,
		  const G4FieldTrack&  CurrentFT, 
                  G4double             requestStep, 
                  G4double             safety,
                  G4int                Step, 
                  G4VPhysicalVolume*   startVolume);
     // Print Method - useful mostly for debugging.

   G4FieldTrack GetEndState();

 public:  // without description

   // void  SetGlobalFieldMgr( G4FieldManager *detectorFieldMgr );
        // The Field Manager of the Detector.

 private:

   G4bool LocateIntersectionPoint( 
	    const  G4FieldTrack&       CurveStartPointTangent,  //  A
	    const  G4FieldTrack&       CurveEndPointTangent,    //  B
	    const  G4ThreeVector&     TrialPoint,              //  E
		   G4FieldTrack&       IntersectPointTangent);  // Output
     // If such an intersection exists, this function 
     // calculate the intersection point of the true path of the particle 
     // with the surface of the current volume (or of one of its daughters). 
     //  (Should use lateral displacement as measure of convergence). 


  //  DATA Members
  // ----------------------------------------------------------------------

 private:

   G4FieldManager *fDetectorFieldMgr; 
     // The  Field Manager of the whole Detector.

   G4Navigator   *fNavigator;

   //  STATE information
   // ------------------

   G4double    fEpsilonStep;
     // Relative accuracy for current Step (Calc.)

   G4FieldTrack    End_PointAndTangent;
     // End point storage

   G4bool      fParticleIsLooping;

   G4int  fVerboseLevel;
     // For debuging purposes

   //  Values for the required accuracies
   //
   G4double  fDelta_One_Step_Value;      //  for one tracking/physics step
   G4double  fDelta_Intersection_Val;    //  for boundary intersection

   //  Their default values ...  (set in G4PropagatemagField.cc)
   //
   static const G4double  fDefault_Delta_One_Step_Value;   // = 0.25 * mm;
   static const G4double  fDefault_Delta_Intersection_Val; // = 0.1 * mm;

   G4int  fmax_loop_count;

   //  Variables to keep track of "abnormal" case - which causes loop
   //
   G4int     fNoZeroStep;                //  Counter of zeroStep
   G4int     fThresholdNo_ZeroSteps;     //  Threshold: above this - action
                       // G4double  fMidPoint_CurveLen_of_LastAttempt= -1;
   G4double  fFull_CurveLen_of_LastAttempt; 
   G4double  fLast_ProposedStepLength; 
};

//  Defines the constructor.
//
#include "G4PropagatorInField.icc"

#endif 
