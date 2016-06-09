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
// $Id: G4PropagatorInField.cc,v 1.42 2008/01/24 08:54:01 gcosmo Exp $
// GEANT4 tag $Name: geant4-09-01-patch-01 $
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
#include "G4GeometryTolerance.hh"
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
  kCarTolerance = G4GeometryTolerance::GetInstance()->GetSurfaceTolerance();

  // In case of too slow progress in finding Intersection Point
  // intermediates Points on the Track must be stored.
  // Initialise the array of Pointers [max_depth+1] to do this  
  
  G4ThreeVector zeroV(0.0,0.0,0.0);
  for (G4int idepth=0; idepth<max_depth+1; idepth++ )
  {
    ptrInterMedFT[ idepth ] = new G4FieldTrack( zeroV, zeroV, 0., 0., 0., 0.);
  }
}

G4PropagatorInField::~G4PropagatorInField()
{
  for ( G4int idepth=0; idepth<max_depth+1; idepth++)
  {
    delete ptrInterMedFT[idepth];
  }
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
  // then return with no action (for now - TODO: some action)
  //
  if(CurrentProposedStepLength<kCarTolerance)
  {
    return kInfinity;
  }

  // Introducing smooth trajectory display (jacek 01/11/2002)
  //
  if (fpTrajectoryFilter)
  {
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
       //G4cout<<"In Locate"<<recalculatedEndPt<<"  and V"<<IntersectPointVelct_G.GetPosition()<<G4endl;
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
     if ( fVerboseLevel > 2 )
       G4cout << " Particle that is stuck will be killed." << G4endl;
     fNoZeroStep = 0; 
  }
  //  G4cout << "G4PropagatorInField returns " << TruePathLength << G4endl;
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
          G4FieldTrack&       IntersectedOrRecalculatedFT, // Out: point found
          G4bool&             recalculatedEndPoint)        // Out: 
{
  // Find Intersection Point ( A, B, E )  of true path AB - start at E.

  G4bool found_approximate_intersection = false;
  G4bool there_is_no_intersection       = false;
  
  G4FieldTrack  CurrentA_PointVelocity = CurveStartPointVelocity; 
  G4FieldTrack  CurrentB_PointVelocity = CurveEndPointVelocity;
  G4ThreeVector CurrentE_Point = TrialPoint;
  G4FieldTrack ApproxIntersecPointV(CurveEndPointVelocity); // FT-Def-Construct
  G4double    NewSafety= -0.0;

  G4bool final_section= true;  // Shows whether current section is last
                               // (i.e. B=full end)
  G4bool first_section=true;
  recalculatedEndPoint= false; 

  G4bool restoredFullEndpoint= false;

  G4int substep_no = 0;
   
  // Limits for substep number
  //
  const G4int max_substeps=   10000;  // Test 120  (old value 100 )
  const G4int warn_substeps=   1000;  //      100  

  // Statistics for substeps
  //
  static G4int max_no_seen= -1; 
  static G4int trigger_substepno_print= warn_substeps - 20 ;

  //--------------------------------------------------------------------------  
  //  Algoritm for the case if progress in founding intersection is too slow.
  //  Process is defined too slow if after N=param_substeps advances on the
  //  path, it will be only 'fraction_done' of the total length.
  //  In this case the remaining length is divided in two half and 
  //  the loop is restarted for each half.  
  //  If progress is still too slow, the division in two halfs continue
  //  until 'max_depth'.
  //--------------------------------------------------------------------------

  const G4int param_substeps=10; // Test value for the maximum number
                                 // of substeps
  const G4double fraction_done=0.3;

  G4bool Second_half=false;      // First half or second half of divided step

  // We need to know this for the 'final_section':
  // real 'final_section' or first half 'final_section'
  // In algorithm it is considered that the 'Second_half' is true
  // and it becomes false only if we are in the first-half of level
  // depthness or if we are in the first section

  G4int depth=0; // Depth counts how many subdivisions of initial step made

#ifdef G4DEBUG_FIELD
  static G4double tolerance= 1.0e-8; 
  G4ThreeVector  StartPosition= CurveStartPointVelocity.GetPosition(); 
  if( (TrialPoint - StartPosition).mag() < tolerance * mm ) 
  {
     G4cerr << "WARNING - G4PropagatorInField::LocateIntersectionPoint()"
            << G4endl
            << "          Intermediate F point is on top of starting point A."
            << G4endl;
     G4Exception("G4PropagatorInField::LocateIntersectionPoint()", 
                 "IntersectionPointIsAtStart", JustWarning,
                 "Intersection point F is exactly at start point A." ); 
  }
#endif

  // Intermediates Points on the Track = Subdivided Points must be stored.
  // Give the initial values to 'InterMedFt'
  // Important is 'ptrInterMedFT[0]', it saves the 'EndCurvePoint'
  //
  *ptrInterMedFT[0] = CurveEndPointVelocity;
  for (G4int idepth=1; idepth<max_depth+1; idepth++ )
  {
    *ptrInterMedFT[idepth]=CurveStartPointVelocity;
  }

  // 'SubStartPoint' is needed to calculate the length of the divided step
  //
  G4FieldTrack SubStart_PointVelocity = CurveStartPointVelocity;
   
  do
  {
    G4int substep_no_p = 0;
    G4bool sub_final_section = false; // the same as final_section,
                                      // but for 'sub_section'
    do // REPEAT param
    {
      G4ThreeVector Point_A = CurrentA_PointVelocity.GetPosition();  
      G4ThreeVector Point_B = CurrentB_PointVelocity.GetPosition();
       
      // F = a point on true AB path close to point E 
      // (the closest if possible)
      //
      ApproxIntersecPointV = GetChordFinder()
                             ->ApproxCurvePointV( CurrentA_PointVelocity, 
                                                  CurrentB_PointVelocity, 
                                                  CurrentE_Point,
                                                  fEpsilonStep );
      //  The above method is the key & most intuitive part ...

#ifdef G4DEBUG_FIELD
      if( ApproxIntersecPointV.GetCurveLength() > 
          CurrentB_PointVelocity.GetCurveLength() * (1.0 + tolerance) )
      {
        G4cerr << "ERROR - G4PropagatorInField::LocateIntersectionPoint()"
               << G4endl
               << "        Intermediate F point is more advanced than"
               << " endpoint B." << G4endl;
        G4Exception("G4PropagatorInField::LocateIntersectionPoint()", 
                    "IntermediatePointConfusion", FatalException,
                    "Intermediate F point is past end B point" ); 
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
        // First relocate to restore any Voxel etc information
        // in the Navigator before calling ComputeStep()
        //
        fNavigator->LocateGlobalPointWithinVolume( Point_A );

        G4ThreeVector PointG;   // Candidate intersection point
        G4double stepLengthAF; 
        G4bool Intersects_AF = IntersectChord( Point_A,   CurrentF_Point,
                                               NewSafety, stepLengthAF,
                                               PointG );
        if( Intersects_AF )
        {
          // G is our new Candidate for the intersection point.
          // It replaces  "E" and we will repeat the test to see if
          // it is a good enough approximate point for us.
          //       B    <- F
          //       E    <- G

          CurrentB_PointVelocity = ApproxIntersecPointV;
          CurrentE_Point = PointG;  

          // By moving point B, must take care if current
          // AF has no intersection to try current FB!!
          //
          final_section= false; 

#ifdef G4VERBOSE
          if( fVerboseLevel > 3 )
          {
            G4cout << "G4PiF::LI> Investigating intermediate point"
                   << " at s=" << ApproxIntersecPointV.GetCurveLength()
                   << " on way to full s="
                   << CurveEndPointVelocity.GetCurveLength() << G4endl;
          }
#endif
        }
        else  // not Intersects_AF
        {  
          // In this case:
          // There is NO intersection of AF with a volume boundary.
          // We must continue the search in the segment FB!
          //
          fNavigator->LocateGlobalPointWithinVolume( CurrentF_Point );

          G4double stepLengthFB;
          G4ThreeVector PointH;

          // Check whether any volumes are encountered by the chord FB
          // ---------------------------------------------------------

          G4bool Intersects_FB = IntersectChord( CurrentF_Point, Point_B, 
                                                 NewSafety, stepLengthFB,
                                                 PointH );
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

            if( final_section  )
            {
              // If B is the original endpoint, this means that whatever
              // volume(s) intersected the original chord, none touch the
              // smaller chords we have used.
              // The value of 'IntersectedOrRecalculatedFT' returned is
              // likely not valid 

              // Check on real final_section or SubEndSection
              //
              if( ((Second_half)&&(depth==0)) || (first_section) )
              {
                there_is_no_intersection = true;   // real final_section
              }
              else
              {
                // end of subsection, not real final section 
                // exit from the and go to the depth-1 level 

                substep_no_p = param_substeps+2;  // exit from the loop

                // but 'Second_half' is still true because we need to find
                // the 'CurrentE_point' for the next loop
                //
                Second_half = true; 
                sub_final_section = true;
           
              }
            }
            else
            {
              // We must restore the original endpoint

              CurrentA_PointVelocity = CurrentB_PointVelocity;  // Got to B
              CurrentB_PointVelocity = CurveEndPointVelocity;
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

        // Change this condition for very strict parameters of propagation 
        //
        if( curveDist*curveDist*(1+2* fEpsilonStep ) < linDistSq )
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

          if( (final_section)&&(Second_half)&&(depth==0) ) // real final section
          {
            recalculatedEndPoint = true;
            IntersectedOrRecalculatedFT = newEndPointFT;
              // So that we can return it, if it is the endpoint!
          }
        }
        if( curveDist < 0.0 )
        {
          G4cerr << "ERROR - G4PropagatorInField::LocateIntersectionPoint()"
                 << G4endl
                 << "        Error in advancing propagation." << G4endl;
          fVerboseLevel = 5; // Print out a maximum of information
          printStatus( CurrentA_PointVelocity,  CurrentB_PointVelocity,
                       -1.0, NewSafety,  substep_no, 0 );
          G4cerr << "        Point A (start) is " << CurrentA_PointVelocity
                 << G4endl;
          G4cerr << "        Point B (end)   is " << CurrentB_PointVelocity
                 << G4endl;
          G4cerr << "        Curve distance is " << curveDist << G4endl;
          G4cerr << G4endl
                 << "The final curve point is not further along"
                 << " than the original!" << G4endl;

          if( recalculatedEndPoint )
          {
            G4cerr << "Recalculation of EndPoint was called with fEpsStep= "
                   << fEpsilonStep << G4endl;
          }
          G4cerr.precision(20);
          G4cerr << " Point A (Curve start)   is " << CurveStartPointVelocity
                 << G4endl;
          G4cerr << " Point B (Curve   end)   is " << CurveEndPointVelocity
                 << G4endl;
          G4cerr << " Point A (Current start) is " << CurrentA_PointVelocity
                 << G4endl;
          G4cerr << " Point B (Current end)   is " << CurrentB_PointVelocity
                 << G4endl;
          G4cerr << " Point S (Sub start)     is " << SubStart_PointVelocity
                 << G4endl;
          G4cerr << " Point E (Trial Point)   is " << CurrentE_Point
                 << G4endl;
          G4cerr << " Point F (Intersection)  is " << ApproxIntersecPointV
                 << G4endl;
          G4cerr << "        LocateIntersection parameters are : Substep no= "
                 << substep_no << G4endl;
          G4cerr << "        Substep depth no= "<< substep_no_p  << " Depth= "
                 << depth << G4endl;

          G4Exception("G4PropagatorInField::LocateIntersectionPoint()",
                      "FatalError", FatalException,
                      "Error in advancing propagation.");
        }

        if(restoredFullEndpoint)
        {
          final_section = restoredFullEndpoint;
          restoredFullEndpoint = false;
        }
      } // EndIf ( E is close enough to the curve, ie point F. )
        // tests ChordAF_Vector.mag() <= maximum_lateral_displacement 

#ifdef G4DEBUG_LOCATE_INTERSECTION  
      if( substep_no >= trigger_substepno_print )
      {
        G4cout << "Difficulty in converging in "
               << "G4PropagatorInField::LocateIntersectionPoint():"
               << G4endl
               << "    Substep no = " << substep_no << G4endl;
        if( substep_no == trigger_substepno_print )
        {
          printStatus( CurveStartPointVelocity, CurveEndPointVelocity,
                       -1.0, NewSafety, 0, 0);
        }
        G4cout << " State of point A: "; 
        printStatus( CurrentA_PointVelocity, CurrentA_PointVelocity,
                     -1.0, NewSafety, substep_no-1, 0);
        G4cout << " State of point B: "; 
        printStatus( CurrentA_PointVelocity, CurrentB_PointVelocity,
                     -1.0, NewSafety, substep_no, 0);
      }
#endif

      substep_no++; 
      substep_no_p++;

    } while (  ( ! found_approximate_intersection )
            && ( ! there_is_no_intersection )     
            && ( substep_no_p <= param_substeps) );  // UNTIL found or
                                                     // failed param substep
    first_section = false;

    if( (!found_approximate_intersection) && (!there_is_no_intersection) )
    {
      G4double did_len = std::abs( CurrentA_PointVelocity.GetCurveLength()
                       - SubStart_PointVelocity.GetCurveLength()); 
      G4double all_len = std::abs( CurrentB_PointVelocity.GetCurveLength()
                       - SubStart_PointVelocity.GetCurveLength());
   
      G4double stepLengthAB;
      G4ThreeVector PointGe;

      // Check if progress is too slow and if it possible to go deeper,
      // then halve the step if so
      //
      if( ( ( did_len )<fraction_done*all_len)
         && (depth<max_depth) && (!sub_final_section) )
      {

        Second_half=false;
        depth++;

        G4double Sub_len = (all_len-did_len)/(2.);
        G4FieldTrack start = CurrentA_PointVelocity;
        G4MagInt_Driver* integrDriver=GetChordFinder()->GetIntegrationDriver();
        integrDriver->AccurateAdvance(start, Sub_len, fEpsilonStep);
        *ptrInterMedFT[depth] = start;
        CurrentB_PointVelocity = *ptrInterMedFT[depth];
 
        // Adjust 'SubStartPoint' to calculate the 'did_length' in next loop
        //
        SubStart_PointVelocity = CurrentA_PointVelocity;

        // Find new trial intersection point needed at start of the loop
        //
        G4ThreeVector Point_A = CurrentA_PointVelocity.GetPosition();
        G4ThreeVector SubE_point = CurrentB_PointVelocity.GetPosition();   
     
        fNavigator->LocateGlobalPointWithinVolume(Point_A);
        G4bool Intersects_AB = IntersectChord(Point_A, SubE_point,
                                              NewSafety, stepLengthAB, PointGe);
        if(Intersects_AB)
        {
          CurrentE_Point = PointGe;
        }
        else
        {
          // No intersection found for first part of curve
          // (CurrentA,InterMedPoint[depth]). Go to the second part
          //
          Second_half = true;
        }
      } // if did_len

      if( (Second_half)&&(depth!=0) )
      {
        // Second part of curve (InterMed[depth],Intermed[depth-1])                       ) 
        // On the depth-1 level normally we are on the 'second_half'

        Second_half = true;

        //  Find new trial intersection point needed at start of the loop
        //
        SubStart_PointVelocity = *ptrInterMedFT[depth];
        CurrentA_PointVelocity = *ptrInterMedFT[depth];
        CurrentB_PointVelocity = *ptrInterMedFT[depth-1];
        G4ThreeVector Point_A    = CurrentA_PointVelocity.GetPosition();
        G4ThreeVector SubE_point = CurrentB_PointVelocity.GetPosition();   
        fNavigator->LocateGlobalPointWithinVolume(Point_A);
        G4bool Intersects_AB = IntersectChord(Point_A, SubE_point, NewSafety,
                                              stepLengthAB, PointGe);
        if(Intersects_AB)
        {
          CurrentE_Point = PointGe;
        }
        else
        {
          final_section = true;
        }
        depth--;
      }
    }  // if(!found_aproximate_intersection)

  } while ( ( ! found_approximate_intersection )
            && ( ! there_is_no_intersection )     
            && ( substep_no <= max_substeps) ); // UNTIL found or failed

  if( substep_no > max_no_seen )
  {
    max_no_seen = substep_no; 
    if( max_no_seen > warn_substeps )
    {
      trigger_substepno_print = max_no_seen-20; // Want to see last 20 steps 
    } 
  }

  if(  ( substep_no >= max_substeps)
      && !there_is_no_intersection
      && !found_approximate_intersection )
  {
    G4cerr << "WARNING - G4PropagatorInField::LocateIntersectionPoint()"
           << G4endl
           << "          Convergence is requiring too many substeps: "
           << substep_no << G4endl;
    G4cerr << "          Abandoning effort to intersect. " << G4endl;
    G4cerr << "          Information on start & current step follows in cout."
           << G4endl;
    G4cout << "WARNING - G4PropagatorInField::LocateIntersectionPoint()"
           << G4endl
           << "          Convergence is requiring too many substeps: "
           << substep_no << G4endl;
    G4cout << "          Found intersection = "
           << found_approximate_intersection << G4endl
           << "          Intersection exists = "
           << !there_is_no_intersection << G4endl;
    G4cout << "          Start and Endpoint of Requested Step:" << G4endl;
    printStatus( CurveStartPointVelocity, CurveEndPointVelocity,
                 -1.0, NewSafety, 0, 0);
    G4cout << G4endl;
    G4cout << "          'Bracketing' starting and endpoint of current Sub-Step"
           << G4endl;
    printStatus( CurrentA_PointVelocity, CurrentA_PointVelocity,
                 -1.0, NewSafety, substep_no-1, 0);
    printStatus( CurrentA_PointVelocity, CurrentB_PointVelocity,
                 -1.0, NewSafety, substep_no, 0);
    G4cout << G4endl;
 
#ifdef FUTURE_CORRECTION
    // Attempt to correct the results of the method // FIX - TODO

    if ( ! found_approximate_intersection )
    { 
      recalculatedEndPoint = true;
      // Return the further valid intersection point -- potentially A ??
      // JA/19 Jan 2006
      IntersectedOrRecalculatedFT = CurrentA_PointVelocity;

      G4cout << "WARNING - G4PropagatorInField::LocateIntersectionPoint()"
             << G4endl
             << "          Did not convergence after " << substep_no
             << " substeps." << G4endl;
      G4cout << "          The endpoint was adjused to pointA resulting"
             << G4endl
             << "          from the last substep: " << CurrentA_PointVelocity
             << G4endl;
    }
#endif

    G4cout.precision( 10 ); 
    G4double done_len = CurrentA_PointVelocity.GetCurveLength(); 
    G4double full_len = CurveEndPointVelocity.GetCurveLength();
    G4cout << "ERROR - G4PropagatorInField::LocateIntersectionPoint()"
           << G4endl
           << "        Undertaken only length: " << done_len
           << " out of " << full_len << " required." << G4endl;
    G4cout << "        Remaining length = " << full_len - done_len << G4endl; 

    G4Exception("G4PropagatorInField::LocateIntersectionPoint()",
                "UnableToLocateIntersection", FatalException,
                "Too many substeps while trying to locate intersection.");
  }
  else if( substep_no >= warn_substeps )
  {  
    int oldprc= G4cout.precision( 10 ); 
    G4cout << "WARNING - G4PropagatorInField::LocateIntersectionPoint()"
           << G4endl
           << "          Undertaken length: "  
           << CurrentB_PointVelocity.GetCurveLength(); 
    G4cout << " - Needed: "  << substep_no << " substeps." << G4endl
           << "          Warning level = " << warn_substeps
           << " and maximum substeps = " << max_substeps << G4endl;
    G4Exception("G4PropagatorInField::LocateIntersectionPoint()",
                "DifficultyToLocateIntersection", JustWarning,
                "Many substeps while trying to locate intersection.");
    G4cout.precision( oldprc ); 
  }
 
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
      || (verboseLevel >= 3) )
  {
    static G4int noPrecision= 4;
    G4cout.precision(noPrecision);
    // G4cout.setf(ios_base::fixed,ios_base::floatfield);
    G4cout << std::setw( 6)  << " " 
           << std::setw( 25) << " Current Position  and  Direction" << " "
           << G4endl; 
    G4cout << std::setw( 5) << "Step#" 
           << std::setw(10) << "  s  " << " "
           << std::setw(10) << "X(mm)" << " "
           << std::setw(10) << "Y(mm)" << " "  
           << std::setw(10) << "Z(mm)" << " "
           << std::setw( 7) << " N_x " << " "
           << std::setw( 7) << " N_y " << " "
           << std::setw( 7) << " N_z " << " " ;
    //            << G4endl; 
    G4cout     // << " >>> "
           << std::setw( 7) << " Delta|N|" << " "
      //   << std::setw( 7) << " Delta(N_z) " << " "
           << std::setw( 9) << "StepLen" << " "  
           << std::setw(12) << "StartSafety" << " "  
           << std::setw( 9) << "PhsStep" << " ";  
    if( startVolume ) {
      G4cout << std::setw(18) << "NextVolume" << " "; 
    }
    G4cout << G4endl;
  }
  if((stepNo == 0) && (verboseLevel <=3)){
     // Recurse to print the start values
     //
     printStatus( StartFT, StartFT, -1.0, safety, -1, startVolume);
   }
   if( verboseLevel <= 3 )
   {
     if( stepNo >= 0)
       G4cout << std::setw( 4) << stepNo << " ";
     else
       G4cout << std::setw( 5) << "Start" ;
     G4cout.precision(8);
     G4cout << std::setw(10) << CurrentFT.GetCurveLength() << " "; 
     G4cout.precision(8);
     G4cout << std::setw(10) << CurrentPosition.x() << " "
            << std::setw(10) << CurrentPosition.y() << " "
            << std::setw(10) << CurrentPosition.z() << " ";
     G4cout.precision(4);
     G4cout << std::setw( 7) << CurrentUnitVelocity.x() << " "
            << std::setw( 7) << CurrentUnitVelocity.y() << " "
            << std::setw( 7) << CurrentUnitVelocity.z() << " ";
     //  G4cout << G4endl; 
     //     G4cout << " >>> " ; 
     G4cout.precision(3); 
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
#if 0
     else
     {
       if( step_len != -1 )
         G4cout << std::setw(12) << "OutOfWorld" << " ";
       else
         G4cout << std::setw(12) << "NotGiven" << " ";
     }
#endif

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

#if 0 
       G4cout << " G4PropagatorInField does not call Navigator::ComputeStep " << G4endl ;
       G4cout << "    step= " << LinearStepLength << " safety= " << NewSafety << G4endl;
       G4cout << "    safety: Origin = " << fPreviousSftOrigin << " val= " << fPreviousSafety << G4endl;
#endif 
    }
    else
    {
       doCallNav= true; 
       // Check whether any volumes are encountered by the chord AB

       // G4cout << " G4PropagatorInField calling Navigator::ComputeStep " << G4endl ;

       LinearStepLength = 
        fNavigator->ComputeStep( StartPointA, ChordAB_Dir,
                                 ChordAB_Length, NewSafety );
       intersects = (LinearStepLength <= ChordAB_Length); 
       // G4Navigator contracts to return k_infinity if len==asked
       // and it did not find a surface boundary at that length
       LinearStepLength = std::min( LinearStepLength, ChordAB_Length);

       // G4cout << " G4PiF got step= " << LinearStepLength << " safety= " << NewSafety << G4endl;

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

    G4cout << " G4PropagatorInField::IntersectChord reports " << G4endl;
    G4cout << " PiF-IC> "
           << "Start="  << std::setw(12) << StartPointA       << " "
           << "End= "   << std::setw(8) << EndPointB         << " "
           << "StepIn=" << std::setw(8) << LinearStepLength  << " "
           << "NewSft=" << std::setw(8) << NewSafety << " " 
           << "CallNav=" << doCallNav      << "  "
           << "Intersects " << intersects     << "  "; 
    if( intersects ) 
      G4cout << "IntrPt=" << std::setw(8) << IntersectionPoint << " " ; 
    G4cout << G4endl;
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
  //                    - CurrentStateA.GetPosition() ).mag2();

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
     if (std::abs(advanceLength)<kCarTolerance)
     {
       advanceLength=(EstimatedEndStateB.GetPosition()
                     -newEndPoint.GetPosition()).mag();
     }
     goodAdvance= 
       integrDriver->AccurateAdvance(newEndPoint, advanceLength, fEpsilonStep);
     //              ***************
  }
  while( !goodAdvance && (++itrial < no_trials) );

  if( goodAdvance )
  {
    retEndPoint= newEndPoint; 
  }
  else
  {
    retEndPoint= EstimatedEndStateB; // Could not improve without major work !!
  }

  //  All the work is done
  //   below are some diagnostics only -- before the return!
  // 
  static const G4String MethodName("G4PropagatorInField::ReEstimateEndpoint");

#ifdef G4VERBOSE
  G4int  latest_good_trials=0;
  if( itrial > 1)
  {
    if( fVerboseLevel > 0 )
    {
      G4cout << MethodName << " called - goodAdv= " << goodAdvance
             << " trials = " << itrial
             << " previous good= " << latest_good_trials
             << G4endl;
    }
    latest_good_trials=0; 
  }
  else
  {   
    latest_good_trials++; 
  }
#endif

#ifdef G4DEBUG_FIELD
  G4double lengthDone = newEndPoint.GetCurveLength() 
                      - CurrentStateA.GetCurveLength(); 
  if( !goodAdvance )
  {
    if( fVerboseLevel >= 3 )
    {
      G4cout << MethodName << "> AccurateAdvance failed " ;
      G4cout << " in " << itrial << " integration trials/steps. " << G4endl;
      G4cout << " It went only " << lengthDone << " instead of " << curveDist
             << " -- a difference of " << curveDist - lengthDone  << G4endl;
      G4cout << " ReEstimateEndpoint> Reset endPoint to original value!"
             << G4endl;
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
  if( goodAdvance )
  {
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
G4PropagatorInField::GimmeTrajectoryVectorAndForgetIt() const
{
  // NB, GimmeThePointsAndForgetThem really forgets them, so it can
  // only be called (exactly) once for each step.

  if (fpTrajectoryFilter)
  {
    return fpTrajectoryFilter->GimmeThePointsAndForgetThem();
  }
  else
  {
    return 0;
  }
}

void 
G4PropagatorInField::SetTrajectoryFilter(G4VCurvedTrajectoryFilter* filter)
{
  fpTrajectoryFilter = filter;
}

void G4PropagatorInField::ClearPropagatorState()
{
  // Goal: Clear all memory of previous steps,  cached information

  fParticleIsLooping= false;
  fNoZeroStep= 0;

  End_PointAndTangent= G4FieldTrack( G4ThreeVector(0.,0.,0.),
                                     G4ThreeVector(0.,0.,0.),
                                     0.0,0.0,0.0,0.0,0.0); 
  fFull_CurveLen_of_LastAttempt = -1; 
  fLast_ProposedStepLength = -1;

  fPreviousSftOrigin= G4ThreeVector(0.,0.,0.);
  fPreviousSafety= 0.0;
}

G4FieldManager* G4PropagatorInField::
FindAndSetFieldManager( G4VPhysicalVolume* pCurrentPhysicalVolume)
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

G4int G4PropagatorInField::SetVerboseLevel( G4int level )
{
  G4int oldval= fVerboseLevel;
  fVerboseLevel= level;

  // Forward the verbose level 'reduced' to ChordFinder,
  // MagIntegratorDriver ... ? 
  //
  G4MagInt_Driver* integrDriver= GetChordFinder()->GetIntegrationDriver(); 
  integrDriver->SetVerboseLevel( fVerboseLevel - 2 ); 
  G4cout << "Set Driver verbosity to " << fVerboseLevel - 2 << G4endl;

  return oldval;
}
