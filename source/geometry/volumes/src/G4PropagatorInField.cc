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
// $Id: G4PropagatorInField.cc,v 1.36 2002-07-29 14:36:35 gcosmo Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
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
// ********************************************************************

#include "G4PropagatorInField.hh"
#include "G4ios.hh"
#include "g4std/iomanip"
 
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
    fEpsilonMinDefault(5.0e-7),
    fEpsilonMaxDefault(0.05),
    fmax_loop_count(10000),
    fNoZeroStep(0)
//    fCharge(0.0), fInitialMomentumModulus(0.0), fMass(0.0)
{
  fEpsilonMin = fEpsilonMinDefault;
  fEpsilonMax = fEpsilonMaxDefault;
  fEpsilonStep = fEpsilonMaxDefault;
  fActionThreshold_NoZeroSteps = 2; 
  fSevereActionThreshold_NoZeroSteps = 10; 
  fAbandonThreshold_NoZeroSteps = 50; 
  fFull_CurveLen_of_LastAttempt = -1; 
  fLast_ProposedStepLength = -1;
  fLargestAcceptableStep = 1000.0 * meter;
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
  // Parameters for adaptive Runge-Kutta integration
  
  G4double      h_TrialStepSize;        // 1st Step Size 
  G4double      TruePathLength = CurrentProposedStepLength;
  G4double      StepTaken = 0.0; 
  G4double      s_length_taken, epsilon ; 
  G4bool        intersects;
  G4bool        first_substep = true;

  G4double      NewSafety;
  fParticleIsLooping = false;

  // Set the field manager if the volume has one, else use the global one
  //
  fCurrentFieldMgr = fDetectorFieldMgr;
  if( pPhysVol)
  {
    G4FieldManager *newFieldMgr = 0;
    newFieldMgr= pPhysVol->GetLogicalVolume()->GetFieldManager(); 
    if ( newFieldMgr ) 
      fCurrentFieldMgr = newFieldMgr;
  }
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
    CurrentProposedStepLength= G4std::min( trialProposedStep,
                                           fLargestAcceptableStep ); 
  }
  epsilon = GetDeltaOneStep() / CurrentProposedStepLength;
  if( epsilon < fEpsilonMin ) epsilon = fEpsilonMin;
  if( epsilon > fEpsilonMax ) epsilon = fEpsilonMax;
  SetEpsilonStep( epsilon );

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
      // WAS = G4Navigator::Locate...

    if( !first_substep)
    {
      fNavigator->LocateGlobalPointWithinVolume( SubStartPoint );
    }

    //   First figure out how far to evolve the particle !
    //   -------------------------------------------------
    //   and (later) with what accuracy to calculate this path.
    //
    h_TrialStepSize = CurrentProposedStepLength - StepTaken;

    //   Next evolve it as far as this allows.
    //   ---------------------------------------
    //
    //   B <- Integrator - limited by the "chord miss" rule.
    //
    s_length_taken = GetChordFinder()->AdvanceChordLimited( 
                             CurrentState,    // Position & velocity
                             h_TrialStepSize,
                             fEpsilonStep);

    fFull_CurveLen_of_LastAttempt = s_length_taken;

    //  On Exit:
    //  CurrentState is updated with the final position and velocity. 
    //
    G4ThreeVector  EndPointB = CurrentState.GetPosition(); 

    // Calculate the direction and length of the chord AB
     
    G4ThreeVector  ChordAB_Vector = EndPointB - SubStartPoint;
    G4double       ChordAB_Length = ChordAB_Vector.mag();  // Magnitude (norm)
    G4ThreeVector  ChordAB_Dir =    ChordAB_Vector.unit();

    // Check whether any volumes are encountered by the chord AB
     
    G4double LinearStepLength = 
        fNavigator->ComputeStep( SubStartPoint, ChordAB_Dir,
                                 ChordAB_Length, NewSafety );
    if( first_substep )
    {
      currentSafety = NewSafety;
    }
    // It might also be possible to update safety in other steps, but
    // it must be done with care.  J.Apostolakis  August 5th, 1997

    intersects = (LinearStepLength <= ChordAB_Length); 
      // G4Navigator contracts to return k_infinity if len==asked
      // and it did not find a surface boundary at that length
    LinearStepLength = G4std::min( LinearStepLength, ChordAB_Length);

    if( intersects )
    {
      // E <- Intersection Point of chord AB and either volume A's surface 
      //                                  or a daughter volume's surface ..

      G4ThreeVector pointE = SubStartPoint + LinearStepLength * ChordAB_Dir;

      G4FieldTrack IntersectPointVelct_G(CurrentState);  // FT-Def-Construct

      // Find the intersection point of AB true path with the surface
      // of vol(A) given our current "estimate" point E. 
  
      G4bool found_intersection = 
        LocateIntersectionPoint( SubStepStartState, CurrentState, 
                                 pointE, IntersectPointVelct_G );
      if( found_intersection )
      {        
        End_PointAndTangent= IntersectPointVelct_G;  // G is our EndPoint ...
        StepTaken = TruePathLength = IntersectPointVelct_G.GetCurveLength()
                                   - OriginalState.GetCurveLength();
         // which is Zero now.
#ifdef G4DEBUG_FIELD
        if( Verbose() > 0 )
          G4cout << " Found intersection after Step of length "
                 << StepTaken << G4endl;
#endif
      }
      else
      {
        intersects= false;  // "Minor" chords do not intersect
      }
    }
    if( !intersects )
    {
      StepTaken += s_length_taken; 
    }
    first_substep = false;

#ifdef G4DEBUG_FIELD
    if( fNoZeroStep > fActionThreshold_NoZeroSteps )
    {
      printStatus( SubStepStartState,  // or OriginalState,
                   CurrentState,
                   CurrentProposedStepLength, 
                   NewSafety,  
                   do_loop_count, 
                   pPhysVol );
    }
#endif

    do_loop_count++;

  } while( (!intersects )
        && (StepTaken + kCarTolerance < CurrentProposedStepLength)  
        && ( do_loop_count < GetMaxLoopCount() ) );

  if( do_loop_count >= GetMaxLoopCount() )
  {
    G4cout << "G4PropagateInField: Warning: Particle is looping - " 
           << " tracking in field will be stopped. " << G4endl;
#ifdef G4VERBOSE
    G4cout << " It has performed " << do_loop_count << " steps in Field " 
           << " while a maximum of " << GetMaxLoopCount() << " are allowed. "
           << G4endl;
#endif
    G4cerr << "G4PropagateInField: Warning: Looping particle. " << G4endl; 
    fParticleIsLooping = true;
  }

  if( !intersects )
  {
    // Chord AB or "minor chords" do not intersect
    // B is the endpoint Step of the current Step.
    // [ But if we were angle limited we could use B 
    //   as a new starting point ? ]
    //  On return we specify the endpoint, point B
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
  if( fabs(OriginalState.GetCurveLength() + TruePathLength 
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
  }
#endif

#ifdef G4DEBUG_FIELD
  if( fNoZeroStep )
  {
    G4cout << " PiF: Step returning= " << StepTaken << G4endl;
    G4cout << " ------------------------------------------------------- "
           << G4endl;
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
#ifdef G4VERBOSE
    G4cout << " WARNING - G4PropagatorInField::ComputeStep():" << G4endl
           << " No progress after "  << fNoZeroStep << " trial steps. "
           << G4endl;
    G4cout << " Particle will be killed." << G4endl; 
#else
    G4cout << " WARNING - G4PropagatorInField::ComputeStep():" << G4endl
           << " Particle that is stuck will be killed." << G4endl;
#endif
    fNoZeroStep = 0; 
    fParticleIsLooping = true;
  }
  return TruePathLength;
}

// --------------------------------------------------------------------------
// G4bool 
// G4PropagatorInField::LocateIntersectionPoint( 
//   const G4FieldTrack&       CurveStartPointVelocity,   //  A
//   const G4FieldTrack&       CurveEndPointVelocity,     //  B
//   const G4ThreeVector&      TrialPoint,                //  E
//         G4FieldTrack&       IntersectPointVelocity)    // Output
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
//    ii) If it returns "false", no intersection was found and 
//       IntersectionPointVelocity is invalid.
// --------------------------------------------------------------------------

G4bool 
G4PropagatorInField::LocateIntersectionPoint( 
  const   G4FieldTrack&       CurveStartPointVelocity,   //  A
  const   G4FieldTrack&       CurveEndPointVelocity,     //  B
  const   G4ThreeVector&      TrialPoint,                //  E
          G4FieldTrack&       IntersectPointVelocity)    // Output
{
  // Find Intersection Point ( A, B, E )  of true path AB - start at E.

  G4bool found_approximate_intersection = false;
  G4bool there_is_no_intersection       = false;

  G4FieldTrack  CurrentA_PointVelocity = CurveStartPointVelocity; 
  G4FieldTrack  CurrentB_PointVelocity = CurveEndPointVelocity;
  G4ThreeVector CurrentE_Point = TrialPoint;

  G4FieldTrack ApproxIntersecPointV(CurveEndPointVelocity); // FT-Def-Construct
  G4double    NewSafety;   
  G4bool      first_step = true;
  G4int       substep_no = 0;

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

    //  The above function is the most difficult part ...
    //        
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
      IntersectPointVelocity = ApproxIntersecPointV;
      IntersectPointVelocity.SetPosition( CurrentE_Point );

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
      // First restore any Voxel etc information in the Navigator
      // before calling ComputeStep   --  not just if (!first_step)
       
      fNavigator->LocateGlobalPointWithinVolume( Point_A );
      // This locate is needed in all cases except for the 
      // original point A, because - presumably - that was 
      // called at the start of the physical Step

      first_step = false;
       
      // Calculate the length and direction of the chord AF
      //
      G4ThreeVector  ChordAF_Vector = CurrentF_Point - Point_A;
      G4double       ChordAF_Length = ChordAF_Vector.mag();  
      G4ThreeVector  ChordAF_Dir    = ChordAF_Vector.unit();

      G4double stepLength = fNavigator->ComputeStep( Point_A, ChordAF_Dir,
                                                     ChordAF_Length, NewSafety);

      G4bool Intersects_AF = (stepLength <= ChordAF_Length);
      stepLength = G4std::min(stepLength, ChordAF_Length);
      if( Intersects_AF )
      {
        // There is an intersection of AF with a volume boundary
        // G <- First Intersection of Chord AF 
        //
        G4ThreeVector PointG = Point_A + stepLength * ChordAF_Dir;

        // G is our new Candidate for the intersection point.
        // It replaces  "E" and we will repeat the test to see if
        // it is a good enough approximate point for us.
        //       B    <- F
        //       E    <- G
        //
        CurrentB_PointVelocity = ApproxIntersecPointV;
        CurrentE_Point = PointG;  
       }
       else  // not Intersects_AF
       {  
         // In this case:
         // There is NO intersection of AF with a volume boundary.
         // We must continue the search in the segment FB!

         // Check whether any volumes are encountered by the chord FB
         // ---------------------------------------------------------
         // Calculate the length and direction of the chord AF
         //
         G4ThreeVector  ChordFB_Vector = Point_B - CurrentF_Point;
         G4double       ChordFB_Length = ChordFB_Vector.mag();  
         G4ThreeVector  ChordFB_Dir    = ChordFB_Vector.unit();

         fNavigator->LocateGlobalPointWithinVolume( CurrentF_Point );
         G4double stepLength =
           fNavigator->ComputeStep( CurrentF_Point, ChordFB_Dir,
                                    ChordFB_Length, NewSafety);

         G4bool Intersects_FB = stepLength <= ChordFB_Length;
         stepLength = G4std::min(stepLength, ChordFB_Length);

         if( Intersects_FB )
         { 
           // There is an intersection of FB with a volume boundary
           // H <- First Intersection of Chord FB 
           //
           G4ThreeVector PointH = CurrentF_Point + stepLength * ChordFB_Dir;

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
           // This means that somehow a volume intersected the original 
           // chord but misses the chord (or series of chords) 
           // we have used.
           // The value of IntersectPointVelocity returned is not valid 
           //
           there_is_no_intersection = true;

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
         G4FieldTrack newEndpoint = CurrentA_PointVelocity;
         GetChordFinder()->GetIntegrationDriver()
           ->AccurateAdvance(newEndpoint, curveDist, fEpsilonStep);
         CurrentB_PointVelocity = newEndpoint;

#ifdef G4DEBUG_FIELD
         static G4int noInaccuracyWarnings = 0; 
         G4int maxNoWarnings = 10;
         if ( (noInaccuracyWarnings < maxNoWarnings ) 
           || (Verbose() > 1) )
         {
           G4cerr << "G4PropagatorInField::LocateIntersectionPoint():"
                  << G4endl
                  << "  Warning: Integration inaccuracy requires " << G4endl
                  << " an adjustment in the step's endpoint "      << G4endl
                  << "   Two mid-points are further apart than their "
                  <<         "curve length difference"             << G4endl 
                  << "   Dist = "       << sqrt(linDistSq)
                  << " curve length = " << curveDist               << G4endl; 
         }
#endif
       }
       if( curveDist < 0.0 )
       {
         G4cerr << "G4PropagatorInField::LocateIntersectionPoint():" << G4endl
                << "Error in advancing propagation." << G4endl;
         printStatus( CurrentA_PointVelocity,  CurrentB_PointVelocity,
                      -1.0, NewSafety,  substep_no, 0);
         G4cout << G4endl
                << "The final curve point is not further along"
                << " than the original!" << G4endl;
         G4Exception("G4PropagatorInField::LocateIntersectionPoint()");
       }
     } // EndIf ( E is close enough to the curve, ie point F. )
       // tests ChordAF_Vector.mag() <= maximum_lateral_displacement 

#ifdef G4VERBOSE
     if( Verbose() > 1 )
        printStatus( CurrentA_PointVelocity,  CurrentB_PointVelocity,
                     -1.0, NewSafety,  substep_no, 0);
#endif

     substep_no++;

   } while (  ( ! found_approximate_intersection )
           && ( ! there_is_no_intersection )     ); // UNTIL found or failed

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
  const G4int verboseLevel=1;
  const G4ThreeVector StartPosition       = StartFT.GetPosition();
  const G4ThreeVector StartUnitVelocity   = StartFT.GetMomentumDir();
  const G4ThreeVector CurrentPosition     = CurrentFT.GetPosition();
  const G4ThreeVector CurrentUnitVelocity = CurrentFT.GetMomentumDir();

  G4double step_len = CurrentFT.GetCurveLength() - StartFT.GetCurveLength();
      
  if( (stepNo == 0) && (verboseLevel <= 3) )
  {
    static G4int noPrecision= 4;
    G4cout.precision(noPrecision);
    // G4cout.setf(ios_base::fixed,ios_base::floatfield);
    G4cout << G4std::setw( 6)  << " " 
           << G4std::setw( 25) << " Current Position  and  Direction" << " "
           << G4endl; 
    G4cout << G4std::setw( 5) << "Step#" << " "
           << G4std::setw( 9) << "X(mm)" << " "
           << G4std::setw( 9) << "Y(mm)" << " "  
           << G4std::setw( 9) << "Z(mm)" << " "
           << G4std::setw( 7) << " N_x " << " "
           << G4std::setw( 7) << " N_y " << " "
           << G4std::setw( 7) << " N_z " << " "
           << G4std::setw( 9) << "StepLen" << " "  
           << G4std::setw(12) << "PhsStep" << " "  
           << G4std::setw(12) << "StartSafety" << " "  
           << G4std::setw(18) << "NextVolume" << " "
           << G4endl;

     // Recurse to print the start values
     //
     printStatus( StartFT, StartFT, -1.0, safety, -1, startVolume);
   }
   if( verboseLevel <= 3 )
   {
     G4cout.precision(3);
     if( stepNo >= 0)
       G4cout << G4std::setw( 5) << stepNo << " ";
     else
       G4cout << G4std::setw( 5) << "Start" << " ";
     G4cout << G4std::setw( 9) << CurrentPosition.x() << " "
            << G4std::setw( 9) << CurrentPosition.y() << " "
            << G4std::setw( 9) << CurrentPosition.z() << " "
            << G4std::setw( 7) << CurrentUnitVelocity.x() << " "
            << G4std::setw( 7) << CurrentUnitVelocity.y() << " "
            << G4std::setw( 7) << CurrentUnitVelocity.z() << " ";
     G4cout << G4std::setw( 9) << step_len << " "; 
     if( requestStep != -1.0 ) 
       G4cout << G4std::setw( 12) << requestStep << " ";
     else
       G4cout << G4std::setw( 12) << "InitialStep" << " "; 
     G4cout << G4std::setw(12) << safety << " ";

     if( startVolume != 0)
     {
       G4cout << G4std::setw(12) << startVolume->GetName() << " ";
     }
     else
     {
       if( step_len != -1 )
         G4cout << G4std::setw(12) << "OutOfWorld" << " ";
       else
         G4cout << G4std::setw(12) << "NotGiven" << " ";
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
                    const G4FieldTrack& aFieldTrack)
{
  G4cout << " PiF: NoZeroStep= " << fNoZeroStep
         << " CurrentProposedStepLength= " << CurrentProposedStepLength
         << " Full_curvelen_last=" << fFull_CurveLen_of_LastAttempt
         << " last proposed step-length= " << fLast_ProposedStepLength 
         << " decreate factor = " << decreaseFactor
         << " step trial = " << stepTrial
         << G4endl;
}
