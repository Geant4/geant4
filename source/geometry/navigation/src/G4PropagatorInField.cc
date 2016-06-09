//
// ********************************************************************
// * DISCLAIMER                                                       *
// *                                                                  *
// * The following disclaimer summarizes all the specific disclaimers *
// * of contributors to this software. The specific disclaimers,which *
// * govern, are listed with their locations in:                      *
// *   http://cern.ch/geant4/license                                  *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.                                                             *
// *                                                                  *
// * This  code  implementation is the  intellectual property  of the *
// * GEANT4 collaboration.                                            *
// * By copying,  distributing  or modifying the Program (or any work *
// * based  on  the Program)  you indicate  your  acceptance of  this *
// * statement, and all its terms.                                    *
// ********************************************************************
//
//
// $Id: G4PropagatorInField.cc,v 1.20 2004/12/02 09:31:23 gcosmo Exp $
// GEANT4 tag $Name: geant4-07-00-cand-03 $
// 
// 
//  This class implements an algorithm to track a particle in a
//  non-uniform magnetic field. It utilises an ODE solver (with
//  the Runge - Kutta method) to evolve the particle, and drives it
//  until the particle has traveled a set distance or it enters a new 
//  volume.
//                                                                     
// 14.10.96 John Apostolakis,   design and implementation
// 17.03.97 John Apostolakis,   renaming new set functions being added
//
// ---------------------------------------------------------------------------

#include "G4PropagatorInField.hh"
#include "G4ios.hh"
#include <iomanip>

#include "G4ThreeVector.hh"
#include "G4VPhysicalVolume.hh"
#include "G4Navigator.hh"
#include "G4VCurvedTrajectoryFilter.hh"
#include "G4ChordFinder.hh"

///////////////////////////////////////////////////////////////////////////
//
// Constructors and destructor

G4PropagatorInField::G4PropagatorInField( G4Navigator    *theNavigator, 
                                          G4FieldManager *detectorFieldMgr )
  : fDetectorFieldMgr(detectorFieldMgr), 
    fCurrentFieldMgr(detectorFieldMgr), 
    fNavigator(theNavigator),
    End_PointAndTangent(G4ThreeVector(0.,0.,0.),
                        G4ThreeVector(0.,0.,0.),0.0,0.0,0.0,0.0,0.0),
    fParticleIsLooping(false),
    fVerboseLevel(0),
    fMax_loop_count(1000),
    fNoZeroStep(0), 
    fCharge(0.0), fInitialMomentumModulus(0.0), fMass(0.0),
    fUseSafetyForOptimisation(true),   // (false) is less sensitive to incorrect safety
    fSetFieldMgr(false),
    fpTrajectoryFilter( 0 )
{
  if(fDetectorFieldMgr) { fEpsilonStep = fDetectorFieldMgr->GetMaximumEpsilonStep();}
  else                  { fEpsilonStep= 1.0e-5; } 
  fActionThreshold_NoZeroSteps = 2; 
  fSevereActionThreshold_NoZeroSteps = 10; 
  fAbandonThreshold_NoZeroSteps = 50; 
  fFull_CurveLen_of_LastAttempt = -1; 
  fLast_ProposedStepLength = -1;
  fLargestAcceptableStep = 1000.0 * meter;

  fPreviousSftOrigin= G4ThreeVector(0.,0.,0.);
  fPreviousSafety= 0.0;
}

G4PropagatorInField::~G4PropagatorInField()
{
}

///////////////////////////////////////////////////////////////////////////
//
// Compute the next geometric Step

G4double
G4PropagatorInField::ComputeStep(
                G4FieldTrack&      pFieldTrack,
                G4double           CurrentProposedStepLength,
                G4double&          currentSafety,                // IN/OUT
                G4VPhysicalVolume* pPhysVol)
{
  // If CurrentProposedStepLength is too small for finding Chords
  // just forget.
  if(CurrentProposedStepLength<kCarTolerance) return DBL_MAX;

  // Introducing smooth trajectory display (jacek 01/11/2002)
  if (fpTrajectoryFilter) {
    fpTrajectoryFilter->CreateNewTrajectorySegment();
  }

  // Parameters for adaptive Runge-Kutta integration
  
  G4double      h_TrialStepSize;        // 1st Step Size 
  G4double      TruePathLength = CurrentProposedStepLength;
  G4double      StepTaken = 0.0; 
  G4double      s_length_taken, epsilon ; 
  G4bool        intersects;
  G4bool        first_substep = true;

  G4double      NewSafety;
  fParticleIsLooping = false;

  // If not yet done, 
  //   Set the field manager to the local  one if the volume has one, 
  //                      or to the global one if not
  //
  if( !fSetFieldMgr ) fCurrentFieldMgr= FindAndSetFieldManager( pPhysVol ); 
  // For the next call, the field manager must again be set
  fSetFieldMgr= false;

  GetChordFinder()->SetChargeMomentumMass(fCharge, fInitialMomentumModulus, fMass);  

  G4FieldTrack  CurrentState(pFieldTrack);
  G4FieldTrack  OriginalState = CurrentState;

  // If the Step length is "infinite", then an approximate-maximum Step
  // length (used to calculate the relative accuracy) must be guessed.
  //
  if( CurrentProposedStepLength >= fLargestAcceptableStep )
  {
    G4ThreeVector StartPointA, VelocityUnit;
    StartPointA  = pFieldTrack.GetPosition();
    VelocityUnit = pFieldTrack.GetMomentumDir();

    G4double trialProposedStep = 1.e2 * ( 10.0 * cm + 
      fNavigator->GetWorldVolume()->GetLogicalVolume()->
                  GetSolid()->DistanceToOut(StartPointA, VelocityUnit) );
    CurrentProposedStepLength= std::min( trialProposedStep,
                                           fLargestAcceptableStep ); 
  }
  epsilon = GetDeltaOneStep() / CurrentProposedStepLength;
  // G4double raw_epsilon= epsilon;
  G4double epsilonMin= fCurrentFieldMgr->GetMinimumEpsilonStep();
  G4double epsilonMax= fCurrentFieldMgr->GetMaximumEpsilonStep();; 
  if( epsilon < epsilonMin ) epsilon = epsilonMin;
  if( epsilon > epsilonMax ) epsilon = epsilonMax;
  SetEpsilonStep( epsilon );

  // G4cout << "G4PiF: Epsilon of current step - raw= " << raw_epsilon
  //        << " final= " << epsilon << G4endl;

  //  Shorten the proposed step in case of earlier problems (zero steps)
  // 
  if( fNoZeroStep > fActionThreshold_NoZeroSteps )
  {
    G4double stepTrial;

    stepTrial= fFull_CurveLen_of_LastAttempt; 
    if( (stepTrial <= 0.0) && (fLast_ProposedStepLength > 0.0) ) 
      stepTrial= fLast_ProposedStepLength; 

    G4double decreaseFactor = 0.9; // Unused default
    if(   (fNoZeroStep < fSevereActionThreshold_NoZeroSteps)
       && (stepTrial > 1000.0*kCarTolerance) )
    {
      // Ensure quicker convergence
      //
      decreaseFactor= 0.1;
    } 
    else
    {
      // We are in significant difficulties, probably at a boundary that
      // is either geometrically sharp or between very different materials.
      // Careful decreases to cope with tolerance are required.
      //
      if( stepTrial > 1000.0*kCarTolerance )
        decreaseFactor = 0.25;     // Try slow decreases
      else if( stepTrial > 100.0*kCarTolerance )
        decreaseFactor= 0.5;       // Try slower decreases
      else if( stepTrial > 10.0*kCarTolerance )
        decreaseFactor= 0.75;      // Try even slower decreases
      else
        decreaseFactor= 0.9;       // Try very slow decreases
     }
     stepTrial *= decreaseFactor;

#ifdef G4DEBUG_FIELD
     PrintStepLengthDiagnostic(CurrentProposedStepLength, decreaseFactor,
                               stepTrial, pFieldTrack);
#endif
     if( stepTrial == 0.0 )
     {
       G4cout << " G4PropagatorInField::ComputeStep "
              << " Particle abandoned due to lack of progress in field."
              << G4endl
              << " Properties : " << pFieldTrack << " "
              << G4endl;
       G4cerr << " G4PropagatorInField::ComputeStep "
              << "  ERROR : attempting a zero step= " << stepTrial << G4endl
              << " while attempting to progress after " << fNoZeroStep
              << " trial steps.  Will abandon step." << G4endl;
         fParticleIsLooping= true;
         return 0;  // = stepTrial;
     }
     if( stepTrial < CurrentProposedStepLength )
       CurrentProposedStepLength = stepTrial;
  }
  fLast_ProposedStepLength = CurrentProposedStepLength;

  G4int do_loop_count = 0; 
  do
  { 
    G4FieldTrack SubStepStartState = CurrentState;
    G4ThreeVector SubStartPoint = CurrentState.GetPosition(); 

    if( !first_substep) {
      fNavigator->LocateGlobalPointWithinVolume( SubStartPoint );
    }

    // How far to attempt to move the particle !
    //
    h_TrialStepSize = CurrentProposedStepLength - StepTaken;

    // Integrate as far as "chord miss" rule allows.
    //
    s_length_taken = GetChordFinder()->AdvanceChordLimited( 
                             CurrentState,    // Position & velocity
                             h_TrialStepSize,
                             fEpsilonStep,
                             fPreviousSftOrigin,
                             fPreviousSafety
                             );
    //  CurrentState is now updated with the final position and velocity. 

    fFull_CurveLen_of_LastAttempt = s_length_taken;

    G4ThreeVector  EndPointB = CurrentState.GetPosition(); 
    G4ThreeVector  InterSectionPointE;
    G4double       LinearStepLength;
 
    // Intersect chord AB with geometry
    intersects= IntersectChord( SubStartPoint, EndPointB, 
				NewSafety,     LinearStepLength, 
				InterSectionPointE );
    // E <- Intersection Point of chord AB and either volume A's surface 
    //                                  or a daughter volume's surface ..

    if( first_substep ) { 
       currentSafety = NewSafety;
    } // Updating safety in other steps is potential future extention

    if( intersects )
    {
       G4FieldTrack IntersectPointVelct_G(CurrentState);  // FT-Def-Construct

       // Find the intersection point of AB true path with the surface
       //   of vol(A), if it exists. Start with point E as first "estimate".
       G4bool recalculatedEndPt= false;
       G4bool found_intersection = 
         LocateIntersectionPoint( SubStepStartState, CurrentState, 
                                  InterSectionPointE, IntersectPointVelct_G,
				  recalculatedEndPt);
       intersects = intersects && found_intersection;
       if( found_intersection ) {        
          End_PointAndTangent= IntersectPointVelct_G;  // G is our EndPoint ...
          StepTaken = TruePathLength = IntersectPointVelct_G.GetCurveLength()
                                      - OriginalState.GetCurveLength();
       } else {
	  // intersects= false;          // "Minor" chords do not intersect
	  if( recalculatedEndPt ){
	     CurrentState= IntersectPointVelct_G; 
	  }
       }
    }
    if( !intersects )
    {
      StepTaken += s_length_taken; 
      // For smooth trajectory display (jacek 01/11/2002)
      if (fpTrajectoryFilter) {
        fpTrajectoryFilter->TakeIntermediatePoint(CurrentState.GetPosition());
      }
    }
    first_substep = false;

#ifdef G4DEBUG_FIELD
    if( fNoZeroStep > fActionThreshold_NoZeroSteps ) {
      printStatus( SubStepStartState,  // or OriginalState,
                   CurrentState,  CurrentProposedStepLength, 
                   NewSafety,     do_loop_count,  pPhysVol );
    }
#endif
#ifdef G4VERBOSE
    if( (fVerboseLevel > 1) && (do_loop_count > fMax_loop_count-10 )) {
      if( do_loop_count == fMax_loop_count-9 ){
	G4cout << "G4PropagatorInField::ComputeStep "
	       << " Difficult track - taking many sub steps." << G4endl;
      }
      printStatus( SubStepStartState, CurrentState, CurrentProposedStepLength, 
		   NewSafety, do_loop_count, pPhysVol );
    }
#endif

    do_loop_count++;

  } while( (!intersects )
        && (StepTaken + kCarTolerance < CurrentProposedStepLength)  
        && ( do_loop_count < fMax_loop_count ) );

  if( do_loop_count >= fMax_loop_count  )
  {
    fParticleIsLooping = true;

    if ( fVerboseLevel > 0 ){
       G4cout << "G4PropagateInField: Killing looping particle " 
              // << " of " << energy  << " energy "
              << " after " << do_loop_count << " field substeps "
              << " totaling " << StepTaken / mm << " mm " ;
       if( pPhysVol )
          G4cout << " in the volume " << pPhysVol->GetName() ; 
       else
         G4cout << " in unknown or null volume. " ; 
       G4cout << G4endl;
    }
  }

  if( !intersects )
  {
    // Chord AB or "minor chords" do not intersect
    // B is the endpoint Step of the current Step.
    //
    End_PointAndTangent = CurrentState; 
    TruePathLength = StepTaken;
  }
  
  // Set pFieldTrack to the return value
  //
  pFieldTrack = End_PointAndTangent;

#ifdef G4VERBOSE
  // Check that "s" is correct
  //
  if( std::fabs(OriginalState.GetCurveLength() + TruePathLength 
      - End_PointAndTangent.GetCurveLength()) > 3.e-4 * TruePathLength )
  {
    G4cerr << " ERROR - G4PropagatorInField::ComputeStep():" << G4endl
           << " Curve length mis-match, is advancement wrong ? " << G4endl;
    G4cerr << " The curve length of the endpoint should be: " 
           << OriginalState.GetCurveLength() + TruePathLength << G4endl
           << " and it is instead: "
           << End_PointAndTangent.GetCurveLength() << "." << G4endl
           << " A difference of: "
           << OriginalState.GetCurveLength() + TruePathLength 
              - End_PointAndTangent.GetCurveLength() << G4endl;
    G4cerr << " Original state= " << OriginalState   << G4endl
	   << " Proposed state= " << End_PointAndTangent << G4endl;
    G4Exception("G4PropagatorInField::ComputeStep()", "IncorrectProposedEndPoint",
		FatalException, 
		"Curve length mis-match between original state and proposed endpoint of propagation.");
  }
#endif

#ifdef G4DEBUG_FIELD
  // static G4std::vector<G4int>  ZeroStepNumberHist(fAbandonThreshold+1);
  if( fNoZeroStep ){
     // ZeroStepNumberHist[fNoZeroStep]++; 
     if( fNoZeroStep > fActionThreshold_NoZeroSteps ){
        G4cout << " PiF: Step returning=" << StepTaken << G4endl;
        G4cout << " ------------------------------------------------------- "
               << G4endl;
     }
  }
#endif

  // In particular anomalous cases, we can get repeated zero steps
  // In order to correct this efficiently, we identify these cases
  // and only take corrective action when they occur.
  // 
  if( TruePathLength < 0.5*kCarTolerance ) 
    fNoZeroStep++;
  else
    fNoZeroStep = 0;

  if( fNoZeroStep > fAbandonThreshold_NoZeroSteps ) { 
     fParticleIsLooping = true;
     G4cout << " WARNING - G4PropagatorInField::ComputeStep():" << G4endl
            << " Zero progress for "  << fNoZeroStep << " attempted steps." 
            << G4endl;
#ifdef G4VERBOSE
     if ( fVerboseLevel > 2 )
       G4cout << " Particle that is stuck will be killed." << G4endl;
#endif
     fNoZeroStep = 0; 
  }

#ifdef G4VERBOSE
  if ( fVerboseLevel > 3 ){
     G4cout << "G4PropagatorInField returns " << TruePathLength << G4endl;
  }
#endif

  return TruePathLength;
}

// --------------------------------------------------------------------------
// G4bool 
// G4PropagatorInField::LocateIntersectionPoint( 
//   const G4FieldTrack&       CurveStartPointVelocity,   //  A
//   const G4FieldTrack&       CurveEndPointVelocity,     //  B
//   const G4ThreeVector&      TrialPoint,                //  E
//         G4FieldTrack&       IntersectedOrRecalculated  // Output
//         G4bool&             recalculated)              // Out
// --------------------------------------------------------------------------
//
// Function that returns the intersection of the true path with the surface
// of the current volume (either the external one or the inner one with one
// of the daughters 
//
//     A = Initial point
//     B = another point 
//
// Both A and B are assumed to be on the true path.
//
//     E is the first point of intersection of the chord AB with 
//     a volume other than A  (on the surface of A or of a daughter)
//
// Convention of Use :
//     i) If it returns "true", then IntersectionPointVelocity is set
//       to the approximate intersection point.
//    ii) If it returns "false", no intersection was found.
//          The validity of IntersectedOrRecalculated depends on 'recalculated'
//        a) if latter is false, then IntersectedOrRecalculated is invalid. 
//        b) if latter is true,  then IntersectedOrRecalculated is
//             the new endpoint, due to a re-integration.
// --------------------------------------------------------------------------

G4bool 
G4PropagatorInField::LocateIntersectionPoint( 
  const   G4FieldTrack&       CurveStartPointVelocity,   //  A
  const   G4FieldTrack&       CurveEndPointVelocity,     //  B
  const   G4ThreeVector&      TrialPoint,                //  E
          G4FieldTrack&       IntersectedOrRecalculatedFT,    // Out: point found
          G4bool&             recalculatedEndPoint)      // Out: 
{
  // Find Intersection Point ( A, B, E )  of true path AB - start at E.

  G4bool found_approximate_intersection = false;
  G4bool there_is_no_intersection       = false;

  G4FieldTrack  CurrentA_PointVelocity = CurveStartPointVelocity; 
  G4FieldTrack  CurrentB_PointVelocity = CurveEndPointVelocity;
  G4ThreeVector CurrentE_Point = TrialPoint;

  G4FieldTrack ApproxIntersecPointV(CurveEndPointVelocity); // FT-Def-Construct
  G4double    NewSafety= -0.0;   
  G4bool final_section= true;  // Shows whether current section is last (ie B=full end)

  recalculatedEndPoint= false; 

  G4bool restoredFullEndpoint= false;

  G4int       substep_no = 0;
  const G4int max_substeps= 100;

  do{                // REPEAT

    G4ThreeVector Point_A = CurrentA_PointVelocity.GetPosition();  
    G4ThreeVector Point_B = CurrentB_PointVelocity.GetPosition();  

    // F = a point on true AB path close to point E  (the closest if possible)
    //
    ApproxIntersecPointV =
      GetChordFinder()->ApproxCurvePointV( CurrentA_PointVelocity, 
                                           CurrentB_PointVelocity, 
                                           CurrentE_Point,
                                           fEpsilonStep );
    //  The above method is the key & most intuitive part ...

#ifdef G4DEBUG_FIELD
    if( ApproxIntersecPointV.GetCurveLength() > 
        CurrentB_PointVelocity.GetCurveLength() * (1.0 + kAngTolerance) ) {
      G4cerr << "Error - Intermediate F point is more advanced than endpoint B." 
	     << G4endl;
      G4Exception("G4PropagatorInField::LocateIntersectionPoint()", 
		  "IntermediatePointConfusion",
		  FatalException, "Intermediate F point is past end B point" ); 
    }
#endif

    G4ThreeVector CurrentF_Point= ApproxIntersecPointV.GetPosition();

    // First check whether EF is small - then F is a good approx. point 
    // Calculate the length and direction of the chord AF
    //
    G4ThreeVector  ChordEF_Vector = CurrentF_Point - CurrentE_Point;

    if ( ChordEF_Vector.mag2() <= sqr(GetDeltaIntersection()) )
    {
      found_approximate_intersection = true;

      // Create the "point" return value
      //
      IntersectedOrRecalculatedFT = ApproxIntersecPointV;
      IntersectedOrRecalculatedFT.SetPosition( CurrentE_Point );

      // Note: in order to return a point on the boundary, 
      //       we must return E. But it is F on the curve.
      //       So we must "cheat": we are using the position at point E
      //       and the velocity at point F !!!
      //
      // This must limit the length we can allow for displacement!

    }
    else  // E is NOT close enough to the curve (ie point F)
    {
      // Check whether any volumes are encountered by the chord AF
      // ---------------------------------------------------------
      // First relocate to restore any Voxel etc information in the Navigator
      //   before calling ComputeStep 
      fNavigator->LocateGlobalPointWithinVolume( Point_A );

      G4ThreeVector PointG;   // Candidate intersection point
      G4double stepLengthAF; 
      G4bool Intersects_AF = IntersectChord( Point_A,   CurrentF_Point,
                                             NewSafety, stepLengthAF,
                                             PointG
                                             );
      if( Intersects_AF )
      {
        // G is our new Candidate for the intersection point.
        // It replaces  "E" and we will repeat the test to see if
        // it is a good enough approximate point for us.
        //       B    <- F
        //       E    <- G
        CurrentB_PointVelocity = ApproxIntersecPointV;
        CurrentE_Point = PointG;  

        // By moving point B, must take care if current AF has no intersection
	//  to try current FB!!
	final_section= false; 

#ifdef G4VERBOSE
	if( fVerboseLevel > 3 ){
	  G4cout << "G4PiF::LI> Investigating intermediate point"
		 << " at s=" << ApproxIntersecPointV.GetCurveLength()
		 << " on way to full s=" << CurveEndPointVelocity.GetCurveLength()
		 << G4endl;
	}
#endif
      }
      else  // not Intersects_AF
      {  
         // In this case:
         // There is NO intersection of AF with a volume boundary.
         // We must continue the search in the segment FB!
         fNavigator->LocateGlobalPointWithinVolume( CurrentF_Point );

         G4double stepLengthFB;
         G4ThreeVector PointH;
         // Check whether any volumes are encountered by the chord FB
         // ---------------------------------------------------------
         G4bool Intersects_FB = 
           IntersectChord( CurrentF_Point, Point_B, 
                           NewSafety,      stepLengthFB,  PointH );
         if( Intersects_FB )
         { 
           // There is an intersection of FB with a volume boundary
           // H <- First Intersection of Chord FB 

           // H is our new Candidate for the intersection point.
           // It replaces  "E" and we will repeat the test to see if
           // it is a good enough approximate point for us.

           // Note that F must be in volume volA  (the same as A)
           // (otherwise AF would meet a volume boundary!)
           //   A    <- F 
           //   E    <- H
           CurrentA_PointVelocity = ApproxIntersecPointV;
           CurrentE_Point = PointH;
         }
         else  // not Intersects_FB
         {
           // There is NO intersection of FB with a volume boundary
	   if( final_section  ){
	     // If B is the original endpoint, this means that whatever volume(s)
	     // intersected the original chord, none touch the smaller chords 
	     // we have used.
	     // The value of IntersectedOrRecalculatedFT returned is likely not valid 
	     //
	     there_is_no_intersection = true;
	   }else{
	     // We must restore the original endpoint
	     CurrentA_PointVelocity= CurrentB_PointVelocity;  // We have got to B
	     CurrentB_PointVelocity= CurveEndPointVelocity;
	     restoredFullEndpoint = true;
	   }

         } // Endif (Intersects_FB)
       } // Endif (Intersects_AF)

       // Ensure that the new endpoints are not further apart in space
       // than on the curve due to different errors in the integration
       //
       G4double linDistSq, curveDist; 
       linDistSq = ( CurrentB_PointVelocity.GetPosition() 
                   - CurrentA_PointVelocity.GetPosition() ).mag2(); 
       curveDist = CurrentB_PointVelocity.GetCurveLength()
                   - CurrentA_PointVelocity.GetCurveLength();
       if( curveDist*(curveDist+2*perMillion ) < linDistSq )
       {
          // Re-integrate to obtain a new B
          //
          G4FieldTrack newEndPointFT=
                  ReEstimateEndpoint( CurrentA_PointVelocity,
                                      CurrentB_PointVelocity,
                                      linDistSq,    // to avoid recalculation
                                      curveDist );
	  G4FieldTrack oldPointVelB = CurrentB_PointVelocity; 
	  CurrentB_PointVelocity = newEndPointFT;

	  if( final_section ){
	     recalculatedEndPoint= true;
	     IntersectedOrRecalculatedFT= newEndPointFT;  // So that we can return it, 
	                                           //  if it is the endpoint!
	  }
       }
       if( curveDist < 0.0 )
       {
         G4cerr << "G4PropagatorInField::LocateIntersectionPoint():" << G4endl
                << "Error in advancing propagation." << G4endl;
	 fVerboseLevel= 5; // Print out a maximum of information
         printStatus( CurrentA_PointVelocity,  CurrentB_PointVelocity,
                      -1.0, NewSafety,  substep_no, 0);
	 G4cerr << " Point A (start) is " << CurrentA_PointVelocity << G4endl;
	 G4cerr << " Point B (end)   is " << CurrentB_PointVelocity << G4endl;
	 G4cerr << " curveDist is " << curveDist << G4endl;
         G4cerr << G4endl
                << "The final curve point is not further along"
                << " than the original!" << G4endl;
         G4Exception("G4PropagatorInField::LocateIntersectionPoint()", "FatalError",
                     FatalException, "Error in advancing propagation.");
       }

       if(restoredFullEndpoint) {
	 final_section= restoredFullEndpoint;	   
	 restoredFullEndpoint=false;
       }

     } // EndIf ( E is close enough to the curve, ie point F. )
       // tests ChordAF_Vector.mag() <= maximum_lateral_displacement 

  } while (  ( ! found_approximate_intersection )
	     && ( ! there_is_no_intersection )     
	     && ( substep_no++ < max_substeps) ); // UNTIL found or failed

#ifdef G4VERBOSE
  if( substep_no >= max_substeps ) {
    G4cerr << "Problem in G4PropagatorInField::LocateIntersectionPoint:"
	   << " Convergence is requiring too many substeps: " << substep_no;
    G4cerr << " Will abandon effort to intersect. " << G4endl;
    G4cerr << " Information on start & current step follows in cout: " << G4endl;
    printStatus( CurrentA_PointVelocity,  CurrentA_PointVelocity,
		 -1.0, NewSafety,  0,          0);
    printStatus( CurrentA_PointVelocity,  CurrentB_PointVelocity,
		 -1.0, NewSafety,  substep_no, 0);
  }
#endif

  return  !there_is_no_intersection; //  Success or failure
}

///////////////////////////////////////////////////////////////////////////
//
// Dumps status of propagator.

void
G4PropagatorInField::printStatus( const G4FieldTrack&        StartFT,
                                  const G4FieldTrack&        CurrentFT, 
                                        G4double             requestStep, 
                                        G4double             safety,
                                        G4int                stepNo, 
                                        G4VPhysicalVolume*   startVolume)
{
  const G4int verboseLevel= fVerboseLevel;
  const G4ThreeVector StartPosition       = StartFT.GetPosition();
  const G4ThreeVector StartUnitVelocity   = StartFT.GetMomentumDir();
  const G4ThreeVector CurrentPosition     = CurrentFT.GetPosition();
  const G4ThreeVector CurrentUnitVelocity = CurrentFT.GetMomentumDir();

  G4double step_len = CurrentFT.GetCurveLength() - StartFT.GetCurveLength();
      
  if( ((stepNo == 0) && (verboseLevel <3))
      || (verboseLevel == 3) )
  {
    static G4int noPrecision= 4;
    G4cout.precision(noPrecision);
    // G4cout.setf(ios_base::fixed,ios_base::floatfield);
    G4cout << std::setw( 6)  << " " 
           << std::setw( 25) << " Current Position  and  Direction" << " "
           << G4endl; 
    G4cout << std::setw( 5) << "Step#" << " "
           << std::setw(10) << "X(mm)" << " "
           << std::setw(10) << "Y(mm)" << " "  
           << std::setw(10) << "Z(mm)" << " "
           << std::setw( 7) << " N_x " << " "
           << std::setw( 7) << " N_y " << " "
           << std::setw( 7) << " N_z " << " "
	   << std::setw( 7) << " Delta|N|" << " "
      //   << std::setw( 7) << " Delta(N_z) " << " "
           << std::setw( 9) << "StepLen" << " "  
           << std::setw(12) << "StartSafety" << " "  
           << std::setw( 9) << "PhsStep" << " "  
           << std::setw(18) << "NextVolume" << " "
           << G4endl;
  }
  if((stepNo == 0) && (verboseLevel <=3)){
     // Recurse to print the start values
     //
     printStatus( StartFT, StartFT, -1.0, safety, -1, startVolume);
   }
   if( verboseLevel <= 3 )
   {
     G4cout.precision(8);
     if( stepNo >= 0)
       G4cout << std::setw( 5) << stepNo << " ";
     else
       G4cout << std::setw( 5) << "Start" << " ";
     G4cout << std::setw(10) << CurrentPosition.x() << " "
            << std::setw(10) << CurrentPosition.y() << " "
            << std::setw(10) << CurrentPosition.z() << " "
            << std::setw( 7) << CurrentUnitVelocity.x() << " "
            << std::setw( 7) << CurrentUnitVelocity.y() << " "
            << std::setw( 7) << CurrentUnitVelocity.z() << " ";
      G4cout.precision(2); 
     G4cout << std::setw( 7) << CurrentFT.GetMomentum().mag()- StartFT.GetMomentum().mag() << " "; 
     //   << std::setw( 7) << CurrentUnitVelocity.z() - InitialUnitVelocity.z() << " ";
     G4cout << std::setw( 9) << step_len << " "; 
     G4cout << std::setw(12) << safety << " ";
     if( requestStep != -1.0 ) 
       G4cout << std::setw( 9) << requestStep << " ";
     else
       G4cout << std::setw( 9) << "Init/NotKnown" << " "; 

     if( startVolume != 0)
     {
       G4cout << std::setw(12) << startVolume->GetName() << " ";
     }
     else
     {
       if( step_len != -1 )
         G4cout << std::setw(12) << "OutOfWorld" << " ";
       else
         G4cout << std::setw(12) << "NotGiven" << " ";
     }

     G4cout << G4endl;
   }
   else // if( verboseLevel > 3 )
   {
     //  Multi-line output
       
     G4cout << "Step taken was " << step_len  
            << " out of PhysicalStep= " <<  requestStep << G4endl;
     G4cout << "Final safety is: " << safety << G4endl;

     G4cout << "Chord length = " << (CurrentPosition-StartPosition).mag()
            << G4endl;
     G4cout << G4endl; 
   }
}

///////////////////////////////////////////////////////////////////////////
//
// Prints Step diagnostics

void 
G4PropagatorInField::PrintStepLengthDiagnostic(
                          G4double CurrentProposedStepLength,
                          G4double decreaseFactor,
                          G4double stepTrial,
                    const G4FieldTrack& )
{
  G4cout << " PiF: NoZeroStep= " << fNoZeroStep
         << " CurrentProposedStepLength= " << CurrentProposedStepLength
         << " Full_curvelen_last=" << fFull_CurveLen_of_LastAttempt
         << " last proposed step-length= " << fLast_ProposedStepLength 
         << " decreate factor = " << decreaseFactor
         << " step trial = " << stepTrial
         << G4endl;
}

G4bool
G4PropagatorInField::IntersectChord( G4ThreeVector  StartPointA, 
                                     G4ThreeVector  EndPointB,
                                     G4double      &NewSafety,
                                     G4double      &LinearStepLength,
                                     G4ThreeVector &IntersectionPoint
                                   )
{
    // Calculate the direction and length of the chord AB
    G4ThreeVector  ChordAB_Vector = EndPointB - StartPointA;
    G4double       ChordAB_Length = ChordAB_Vector.mag();  // Magnitude (norm)
    G4ThreeVector  ChordAB_Dir =    ChordAB_Vector.unit();
    G4bool intersects;

    G4ThreeVector OriginShift = StartPointA - fPreviousSftOrigin ;
    G4double      MagSqShift  = OriginShift.mag2() ;
    G4double      currentSafety;
    G4bool        doCallNav= false;

    if( MagSqShift >= sqr(fPreviousSafety) )
    {
        currentSafety = 0.0 ;
    }else{
        currentSafety = fPreviousSafety - std::sqrt(MagSqShift) ;
    }

    if( fUseSafetyForOptimisation && (ChordAB_Length <= currentSafety) )
    {
       // The Step is guaranteed to be taken

       LinearStepLength = ChordAB_Length;
       intersects = false;

       NewSafety= currentSafety;
    }
    else
    {
       doCallNav= true; 
       // Check whether any volumes are encountered by the chord AB
       LinearStepLength = 
        fNavigator->ComputeStep( StartPointA, ChordAB_Dir,
                                 ChordAB_Length, NewSafety );
       intersects = (LinearStepLength <= ChordAB_Length); 
       // G4Navigator contracts to return k_infinity if len==asked
       // and it did not find a surface boundary at that length
       LinearStepLength = std::min( LinearStepLength, ChordAB_Length);

       // Save the last calculated safety!
       fPreviousSftOrigin = StartPointA;
       fPreviousSafety= NewSafety;

       if( intersects ){
          // Intersection Point of chord AB and either volume A's surface 
          //                                or a daughter volume's surface ..
          IntersectionPoint = StartPointA + LinearStepLength * ChordAB_Dir;
       }
    }

#ifdef DEBUG_INTERSECTS_CHORD
    // printIntersection( 
    // StartPointA, EndPointB, LinearStepLength, IntersectionPoint, NewSafety

    G4cout << "Start="  << std::setw(12) << StartPointA       << " "
           << "End= "   << std::setw(8) << EndPointB         << " "
           << "StepIn=" << std::setw(8) << LinearStepLength  << " "
           << "NewSft=" << std::setw(8) << NewSafety
           << "NavCall" << doCallNav      << "  "
           << "In T/F " << intersects     << "  " 
           << "IntrPt=" << std::setw(8) << IntersectionPoint << " " 
           << G4endl;
#endif

    return intersects;
}

// --------------------- oooo000000000000oooo ----------------------------

G4FieldTrack G4PropagatorInField::
ReEstimateEndpoint( const G4FieldTrack &CurrentStateA,  
                    const G4FieldTrack &EstimatedEndStateB,
                          G4double      linearDistSq,
                          G4double      curveDist
                  )
{
  // G4double checkCurveDist= EstimatedEndStateB.GetCurveLength()
  //   - CurrentStateA.GetCurveLength();
  // G4double checkLinDistSq= (EstimatedEndStateB.GetPosition()
  //		    - CurrentStateA.GetPosition() ).mag2();

  G4FieldTrack newEndPoint( CurrentStateA );
  G4MagInt_Driver* integrDriver= GetChordFinder()->GetIntegrationDriver(); 

  G4FieldTrack retEndPoint( CurrentStateA );
  G4bool goodAdvance;
  G4int  itrial=0;
  const G4int no_trials= 20;

  G4double endCurveLen= EstimatedEndStateB.GetCurveLength();
  do
  {
     G4double currentCurveLen= newEndPoint.GetCurveLength();
     G4double advanceLength= endCurveLen - currentCurveLen ; 

     goodAdvance= 
       integrDriver->AccurateAdvance(newEndPoint, advanceLength, fEpsilonStep);
     //              ***************
  }
  while( !goodAdvance && (++itrial < no_trials) );

  if( goodAdvance ) {
    retEndPoint= newEndPoint; 
  }else{
    retEndPoint= EstimatedEndStateB; // Could not improve without major work !!
  }

  //  All the work is done
  //   below are some diagnostics only -- before the return!
  // 
  static const G4String MethodName("G4PropagatorInField::ReEstimateEndpoint");
#ifdef G4VERBOSE
  G4int  latest_good_trials=0;
  if( itrial > 1) {
    if( fVerboseLevel > 0 ) {
      G4cout << MethodName << " called - goodAdv= " << goodAdvance
	 << " trials = " << itrial << " previous good= " << latest_good_trials
	 << G4endl;
    }
    latest_good_trials=0; 
  }else{   
    latest_good_trials++; 
  }
#endif

#ifdef G4DEBUG_FIELD
  G4double lengthDone=  newEndPoint.GetCurveLength() 
                           - CurrentStateA.GetCurveLength(); 
  if( !goodAdvance ) {
    if( fVerboseLevel >= 3 ){
      G4cout << MethodName << "> AccurateAdvance failed " ;
      G4cout << " in " << itrial << " integration trials/steps. " << G4endl
      G4cout << " It went only " << lengthDone << " instead of " << curveDist 
	     << " -- a difference of " << curveDist - lengthDone  << G4endl;
      G4cout << " ReEstimateEndpoint> Reset endPoint to original value!" << G4endl;
    }
  }

  static G4int noInaccuracyWarnings = 0; 
  G4int maxNoWarnings = 10;
  if (  (noInaccuracyWarnings < maxNoWarnings ) 
       || (fVerboseLevel > 1) )
    {
      G4cerr << "G4PropagatorInField::LocateIntersectionPoint():"
             << G4endl
             << " Warning: Integration inaccuracy requires" 
             <<   " an adjustment in the step's endpoint."  << G4endl
             << "   Two mid-points are further apart than their"
             <<   " curve length difference"                << G4endl 
             << "   Dist = "       << std::sqrt(linearDistSq)
             << " curve length = " << curveDist             << G4endl; 
      G4cerr << " Correction applied is " 
             << (newEndPoint.GetPosition()-EstimatedEndStateB.GetPosition()).mag()
             << G4endl;
    }
#else
  // Statistics on the RMS value of the corrections
  static G4int    noCorrections=0;
  static G4double sumCorrectionsSq = 0;
  noCorrections++; 
  if( goodAdvance ){
    sumCorrectionsSq += (EstimatedEndStateB.GetPosition() - 
			 newEndPoint.GetPosition()).mag2();
  }
  linearDistSq -= curveDist; // To use linearDistSq ... !
#endif

  return retEndPoint;
}

// Access the points which have passed through the filter. The
// points are stored as ThreeVectors for the initial impelmentation
// only (jacek 30/10/2002)
// Responsibility for deleting the points lies with
// SmoothTrajectoryPoint, which is the points' final
// destination. The points pointer is set to NULL, to ensure that
// the points are not re-used in subsequent steps, therefore THIS
// METHOD MUST BE CALLED EXACTLY ONCE PER STEP. (jacek 08/11/2002)


std::vector<G4ThreeVector>*
G4PropagatorInField::GimmeTrajectoryVectorAndForgetIt() const {
  // NB, GimmeThePointsAndForgetThem really forgets them, so it can
  // only be called (exactly) once for each step.
  if (fpTrajectoryFilter) {
    return fpTrajectoryFilter->GimmeThePointsAndForgetThem();
  } else {
    return NULL;
  }
}

void 
G4PropagatorInField::SetTrajectoryFilter(G4VCurvedTrajectoryFilter* filter) {
  fpTrajectoryFilter = filter;
}


void G4PropagatorInField::ClearPropagatorState()
{
  G4Exception("G4PropagatorInField::ClearPropagatorState()", "NotImplemented",
              FatalException, "Functionality not yet implemented.");
}

G4FieldManager* 
G4PropagatorInField::FindAndSetFieldManager( G4VPhysicalVolume* pCurrentPhysicalVolume)
{
  G4FieldManager* currentFieldMgr;

  currentFieldMgr = fDetectorFieldMgr;
  if( pCurrentPhysicalVolume)
  {
     G4FieldManager *newFieldMgr = 0;
     newFieldMgr= pCurrentPhysicalVolume->GetLogicalVolume()->GetFieldManager();
     if ( newFieldMgr ) 
        currentFieldMgr = newFieldMgr;
  }
  fCurrentFieldMgr= currentFieldMgr;

  // Flag that field manager has been set.
  fSetFieldMgr= true;

  return currentFieldMgr;
}

G4int G4PropagatorInField::SetVerboseLevel( G4int Verbose )
{
  G4int oldval= fVerboseLevel;
  fVerboseLevel= Verbose;

  // Forward the verbose level 'reduced' to ChordFinder, MagIntegratorDriver ... ? 
  G4MagInt_Driver* integrDriver= GetChordFinder()->GetIntegrationDriver(); 
  integrDriver->SetVerboseLevel( Verbose - 2 ); 
  G4cout << "Set Driver verbosity to " << Verbose - 2 << G4endl;

  return oldval;
}
