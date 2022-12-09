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
// class G4PropagatorInField Implementation
// 
//  This class implements an algorithm to track a particle in a
//  non-uniform magnetic field. It utilises an ODE solver (with
//  the Runge - Kutta method) to evolve the particle, and drives it
//  until the particle has traveled a set distance or it enters a new 
//  volume.
//                                                                     
// 14.10.96 John Apostolakis, design and implementation
// 17.03.97 John Apostolakis, renaming new set functions being added
// ---------------------------------------------------------------------------

#include <iomanip>

#include "G4PropagatorInField.hh"
#include "G4ios.hh"
#include "G4SystemOfUnits.hh"
#include "G4ThreeVector.hh"
#include "G4Material.hh"
#include "G4VPhysicalVolume.hh"
#include "G4Navigator.hh"
#include "G4GeometryTolerance.hh"
#include "G4VCurvedTrajectoryFilter.hh"
#include "G4ChordFinder.hh"
#include "G4MultiLevelLocator.hh"

// ---------------------------------------------------------------------------
// Constructors and destructor
//
G4PropagatorInField::G4PropagatorInField( G4Navigator* theNavigator, 
                                          G4FieldManager* detectorFieldMgr,
                                          G4VIntersectionLocator* vLocator  )
  : fDetectorFieldMgr(detectorFieldMgr), 
    fNavigator(theNavigator),
    fCurrentFieldMgr(detectorFieldMgr),
    End_PointAndTangent(G4ThreeVector(0.,0.,0.),
                        G4ThreeVector(0.,0.,0.),0.0,0.0,0.0,0.0,0.0)
{
  fEpsilonStep = (fDetectorFieldMgr != nullptr)
               ? fDetectorFieldMgr->GetMaximumEpsilonStep() : 1.0e-5;
  fLargestAcceptableStep = 1000.0 * meter;

  fPreviousSftOrigin = G4ThreeVector(0.,0.,0.);
  kCarTolerance = G4GeometryTolerance::GetInstance()->GetSurfaceTolerance();
  fZeroStepThreshold = std::max( 1.0e5 * kCarTolerance, 1.0e-1 * micrometer );

#ifdef G4DEBUG_FIELD
  G4cout << " PiF: Zero Step Threshold set to "
         << fZeroStepThreshold / millimeter
         << " mm." << G4endl;
  G4cout << " PiF:   Value of kCarTolerance = "
         << kCarTolerance / millimeter 
         << " mm. " << G4endl;
  fVerboseLevel = 2;
  fVerbTracePiF = true;   
#endif 

  // Defining Intersection Locator and his parameters
  if ( vLocator == nullptr )
  {
    fIntersectionLocator = new G4MultiLevelLocator(theNavigator);
    fAllocatedLocator = true;
  }
  else
  {
    fIntersectionLocator = vLocator;
    fAllocatedLocator = false;
  }
  RefreshIntersectionLocator();  //  Copy all relevant parameters
}

// ---------------------------------------------------------------------------
//
G4PropagatorInField::~G4PropagatorInField()
{
  if(fAllocatedLocator)  { delete  fIntersectionLocator; }
}

// ---------------------------------------------------------------------------
// Update the IntersectionLocator with current parameters
//
void G4PropagatorInField::RefreshIntersectionLocator()
{
  fIntersectionLocator->SetEpsilonStepFor(fEpsilonStep);
  fIntersectionLocator->SetDeltaIntersectionFor(fCurrentFieldMgr->GetDeltaIntersection());
  fIntersectionLocator->SetChordFinderFor(GetChordFinder());
  fIntersectionLocator->SetSafetyParametersFor( fUseSafetyForOptimisation);
}

// ---------------------------------------------------------------------------
// Compute the next geometric Step
//
G4double G4PropagatorInField::ComputeStep(
                G4FieldTrack&      pFieldTrack,
                G4double           CurrentProposedStepLength,
                G4double&          currentSafety,                // IN/OUT
                G4VPhysicalVolume* pPhysVol,
                G4bool             canRelaxDeltaChord)
{  
  GetChordFinder()->OnComputeStep();
  const G4double deltaChord = GetChordFinder()->GetDeltaChord();

  // If CurrentProposedStepLength is too small for finding Chords
  // then return with no action (for now - TODO: some action)
  //
  const char* methodName = "G4PropagatorInField::ComputeStep";
  if (CurrentProposedStepLength<kCarTolerance)
  {
    return kInfinity;
  }

  // Introducing smooth trajectory display (jacek 01/11/2002)
  //
  if (fpTrajectoryFilter)
  {
    fpTrajectoryFilter->CreateNewTrajectorySegment();
  }

  fFirstStepInVolume = fNewTrack ? true : fLastStepInVolume;
  fLastStepInVolume = false;
  fNewTrack = false; 

  if( fVerboseLevel > 2 )
  {
    G4cout << methodName << " called" << G4endl;
    G4cout << "   Starting FT: " << pFieldTrack;
    G4cout << "   Requested length = " << CurrentProposedStepLength << G4endl;
    G4cout << "   PhysVol = ";
    if( pPhysVol )
       G4cout << pPhysVol->GetName() << G4endl;
    else
       G4cout << " N/A ";
    G4cout << G4endl;
  }
  
  // Parameters for adaptive Runge-Kutta integration
  
  G4double h_TrialStepSize;        // 1st Step Size 
  G4double TruePathLength = CurrentProposedStepLength;
  G4double StepTaken = 0.0; 
  G4double s_length_taken, epsilon; 
  G4bool   intersects;
  G4bool   first_substep = true;

  G4double NewSafety;
  fParticleIsLooping = false;

  // If not yet done, 
  //   Set the field manager to the local  one if the volume has one, 
  //                      or to the global one if not
  //
  if( !fSetFieldMgr )
  {
    fCurrentFieldMgr = FindAndSetFieldManager( pPhysVol );
  }
  fSetFieldMgr = false; // For next call, the field manager must be set again

  G4FieldTrack CurrentState(pFieldTrack);
  G4FieldTrack OriginalState = CurrentState;

  // If the Step length is "infinite", then an approximate-maximum Step
  // length (used to calculate the relative accuracy) must be guessed
  //
  if( CurrentProposedStepLength >= fLargestAcceptableStep )
  {
    G4ThreeVector StartPointA, VelocityUnit;
    StartPointA  = pFieldTrack.GetPosition();
    VelocityUnit = pFieldTrack.GetMomentumDir();

    G4double trialProposedStep = 1.e2 * ( 10.0 * cm + 
      fNavigator->GetWorldVolume()->GetLogicalVolume()->
                  GetSolid()->DistanceToOut(StartPointA, VelocityUnit) );
    CurrentProposedStepLength = std::min( trialProposedStep,
                                          fLargestAcceptableStep ); 
  }
  epsilon = fCurrentFieldMgr->GetDeltaOneStep() / CurrentProposedStepLength;
  G4double epsilonMin= fCurrentFieldMgr->GetMinimumEpsilonStep();
  G4double epsilonMax= fCurrentFieldMgr->GetMaximumEpsilonStep();
  if( epsilon < epsilonMin )  { epsilon = epsilonMin; }
  if( epsilon > epsilonMax )  { epsilon = epsilonMax; }
  SetEpsilonStep( epsilon );

  // Values for Intersection Locator has to be updated on each call for the
  // case that CurrentFieldManager has changed from the one of previous step
  //
  RefreshIntersectionLocator();

  // Shorten the proposed step in case of earlier problems (zero steps)
  // 
  if( fNoZeroStep > fActionThreshold_NoZeroSteps )
  {
    G4double stepTrial;

    stepTrial = fFull_CurveLen_of_LastAttempt; 
    if( (stepTrial <= 0.0) && (fLast_ProposedStepLength > 0.0) )
    {
      stepTrial = fLast_ProposedStepLength; 
    }

    G4double decreaseFactor = 0.9; // Unused default
    if(   (fNoZeroStep < fSevereActionThreshold_NoZeroSteps)
       && (stepTrial > 100.0*fZeroStepThreshold) )
    {
      // Attempt quick convergence
      //
      decreaseFactor= 0.25;
    } 
    else
    {
      // We are in significant difficulties, probably at a boundary that
      // is either geometrically sharp or between very different materials.
      // Careful decreases to cope with tolerance are required
      //
      if( stepTrial > 100.0*fZeroStepThreshold )
        decreaseFactor = 0.35;     // Try decreasing slower
      else if( stepTrial > 30.0*fZeroStepThreshold )
        decreaseFactor= 0.5;       // Try yet slower decrease
      else if( stepTrial > 10.0*fZeroStepThreshold )
        decreaseFactor= 0.75;      // Try even slower decreases
      else
        decreaseFactor= 0.9;       // Try very slow decreases
     }
     stepTrial *= decreaseFactor;

#ifdef G4DEBUG_FIELD
     if( fVerboseLevel > 2
      || (fNoZeroStep >= fSevereActionThreshold_NoZeroSteps) )
     {
        G4cerr << " " << methodName
               << "  Decreasing step after " << fNoZeroStep << " zero steps "
               << " - in volume " << pPhysVol;
        if( pPhysVol )
           G4cerr << " with name " << pPhysVol->GetName();
        else
           G4cerr << " i.e. *unknown* volume.";
        G4cerr << G4endl;
        PrintStepLengthDiagnostic(CurrentProposedStepLength, decreaseFactor,
                                  stepTrial, pFieldTrack);
     }
#endif
     if( stepTrial == 0.0 )  //  Change to make it < 0.1 * kCarTolerance ??
     {
       std::ostringstream message;
       message << "Particle abandoned due to lack of progress in field."
               << G4endl
               << "  Properties : " << pFieldTrack << G4endl
               << "  Attempting a zero step = " << stepTrial << G4endl
               << "  while attempting to progress after " << fNoZeroStep
               << " trial steps. Will abandon step.";
       G4Exception(methodName, "GeomNav1002", JustWarning, message);
       fParticleIsLooping = true;
       return 0;  // = stepTrial;
     }
     if( stepTrial < CurrentProposedStepLength )
     {
       CurrentProposedStepLength = stepTrial;
     }
  }
  fLast_ProposedStepLength = CurrentProposedStepLength;

  G4int do_loop_count = 0; 
  do  // Loop checking, 07.10.2016, JA
  { 
    G4FieldTrack SubStepStartState = CurrentState;
    G4ThreeVector SubStartPoint = CurrentState.GetPosition(); 
    
    if(!first_substep)
    {
      if( fVerboseLevel > 4 )
      {
        G4cout << " PiF: Calling Nav/Locate Global Point within-Volume "
               << G4endl;
      }
      fNavigator->LocateGlobalPointWithinVolume( SubStartPoint );
    }

    // How far to attempt to move the particle !
    //
    h_TrialStepSize = CurrentProposedStepLength - StepTaken;

    if (canRelaxDeltaChord &&
        fIncreaseChordDistanceThreshold > 0  &&
        do_loop_count > fIncreaseChordDistanceThreshold && 
        do_loop_count % fIncreaseChordDistanceThreshold == 0)
    {
        GetChordFinder()->SetDeltaChord(
          GetChordFinder()->GetDeltaChord() * 2.0
        );
    }

    // Integrate as far as "chord miss" rule allows.
    //
    s_length_taken = GetChordFinder()->AdvanceChordLimited( 
                             CurrentState,    // Position & velocity
                             h_TrialStepSize,
                             fEpsilonStep,
                             fPreviousSftOrigin,
                             fPreviousSafety );
      // CurrentState is now updated with the final position and velocity

    fFull_CurveLen_of_LastAttempt = s_length_taken;

    G4ThreeVector EndPointB = CurrentState.GetPosition(); 
    G4ThreeVector InterSectionPointE;
    G4double      LinearStepLength;
 
    // Intersect chord AB with geometry
    //
    intersects= IntersectChord( SubStartPoint, EndPointB, 
                                NewSafety, LinearStepLength, 
                                InterSectionPointE );
      // E <- Intersection Point of chord AB and either volume A's surface 
      //                                  or a daughter volume's surface ..

    if( first_substep )
    { 
       currentSafety = NewSafety;
    } // Updating safety in other steps is potential future extention

    if( intersects )
    {
       G4FieldTrack IntersectPointVelct_G(CurrentState);  // FT-Def-Construct

       // Find the intersection point of AB true path with the surface
       //   of vol(A), if it exists. Start with point E as first "estimate".
       G4bool recalculatedEndPt = false;
       
       G4bool found_intersection = fIntersectionLocator->
         EstimateIntersectionPoint( SubStepStartState, CurrentState, 
                                    InterSectionPointE, IntersectPointVelct_G,
                                    recalculatedEndPt, fPreviousSafety,
                                    fPreviousSftOrigin);
       intersects = found_intersection;
       if( found_intersection )
       {        
          End_PointAndTangent= IntersectPointVelct_G;  // G is our EndPoint ...
          StepTaken = TruePathLength = IntersectPointVelct_G.GetCurveLength()
                                     - OriginalState.GetCurveLength();
       }
       else
       {
          // Either "minor" chords do not intersect
          // or else stopped (due to too many steps)
          //
          if( recalculatedEndPt )
          {
             G4double endAchieved = IntersectPointVelct_G.GetCurveLength();
             G4double endExpected = CurrentState.GetCurveLength(); 

             // Detect failure - due to too many steps
             G4bool shortEnd = endAchieved
                             < (endExpected*(1.0-CLHEP::perMillion));

             G4double stepAchieved = endAchieved
                                   - SubStepStartState.GetCurveLength();

             // Update remaining state - must work for 'full' step or
             // abandonned intersection
             //
             CurrentState = IntersectPointVelct_G;
             s_length_taken = stepAchieved;
             if( shortEnd )
             {
                fParticleIsLooping = true;
             } 
          }
       }
    }
    if( !intersects )
    {
      StepTaken += s_length_taken;

      if (fpTrajectoryFilter) // For smooth trajectory display (jacek 1/11/2002)
      {
        fpTrajectoryFilter->TakeIntermediatePoint(CurrentState.GetPosition());
      }
    }
    first_substep = false;

#ifdef G4DEBUG_FIELD
    if( fNoZeroStep > fActionThreshold_NoZeroSteps )
    {
      if( fNoZeroStep > fSevereActionThreshold_NoZeroSteps )
        G4cout << " Above 'Severe Action' threshold -- for Zero steps.  ";
      else
        G4cout << " Above 'action' threshold -- for Zero steps.  ";         
      G4cout << " Number of zero steps = " << fNoZeroStep << G4endl;
      printStatus( SubStepStartState,  // or OriginalState,
                   CurrentState, CurrentProposedStepLength, 
                   NewSafety, do_loop_count, pPhysVol );
    }
    if( (fVerboseLevel > 1) && (do_loop_count > fMax_loop_count-10 ))
    {
      if( do_loop_count == fMax_loop_count-9 )
      {
        G4cout << " G4PropagatorInField::ComputeStep(): " << G4endl
               << "  Difficult track - taking many sub steps." << G4endl;
        printStatus( SubStepStartState, SubStepStartState, CurrentProposedStepLength, 
                     NewSafety, 0, pPhysVol );        
      }
      printStatus( SubStepStartState, CurrentState, CurrentProposedStepLength, 
                   NewSafety, do_loop_count, pPhysVol );
    }
#endif

    ++do_loop_count;

  } while( (!intersects )
        && (!fParticleIsLooping)
        && (StepTaken + kCarTolerance < CurrentProposedStepLength)  
        && ( do_loop_count < fMax_loop_count ) );

  if(  do_loop_count >= fMax_loop_count
    && (StepTaken + kCarTolerance < CurrentProposedStepLength) )
  {
    fParticleIsLooping = true;
  }
  if ( ( fParticleIsLooping ) && (fVerboseLevel > 0) )
  {
    ReportLoopingParticle( do_loop_count, StepTaken,
                           CurrentProposedStepLength, methodName,
                           CurrentState.GetMomentum(), pPhysVol );
  }
    
  if( !intersects )
  {
    // Chord AB or "minor chords" do not intersect
    // B is the endpoint Step of the current Step.
    //
    End_PointAndTangent = CurrentState; 
    TruePathLength = StepTaken;   //  Original code
 
    // Tried the following to avoid potential issue with round-off error
    // - but has issues... Suppressing this change JA 2015/05/02
    // TruePathLength = CurrentProposedStepLength;
  }
  fLastStepInVolume = intersects;
  
  // Set pFieldTrack to the return value
  //
  pFieldTrack = End_PointAndTangent;

#ifdef G4VERBOSE
  // Check that "s" is correct
  //
  if( std::fabs(OriginalState.GetCurveLength() + TruePathLength 
      - End_PointAndTangent.GetCurveLength()) > 3.e-4 * TruePathLength )
  {
    std::ostringstream message;
    message << "Curve length mis-match between original state "
            << "and proposed endpoint of propagation." << G4endl
            << "  The curve length of the endpoint should be: " 
            << OriginalState.GetCurveLength() + TruePathLength << G4endl
            << "  and it is instead: "
            << End_PointAndTangent.GetCurveLength() << "." << G4endl
            << "  A difference of: "
            << OriginalState.GetCurveLength() + TruePathLength 
               - End_PointAndTangent.GetCurveLength() << G4endl
            << "  Original state = " << OriginalState   << G4endl
            << "  Proposed state = " << End_PointAndTangent;
    G4Exception(methodName, "GeomNav0003", FatalException, message);
  }
#endif

  if( TruePathLength+kCarTolerance >= CurrentProposedStepLength )
  {
     fNoZeroStep = 0;     
  }
  else
  {     
     // In particular anomalous cases, we can get repeated zero steps
     // We identify these cases and take corrective action when they occur.
     // 
     if( TruePathLength < std::max( fZeroStepThreshold, 0.5*kCarTolerance ) )
     {
        ++fNoZeroStep;
     }
     else
     {
        fNoZeroStep = 0;
     }
  }
  if( fNoZeroStep > fAbandonThreshold_NoZeroSteps )
  { 
     fParticleIsLooping = true;
     ReportStuckParticle( fNoZeroStep, CurrentProposedStepLength,
                          fFull_CurveLen_of_LastAttempt, pPhysVol );
     fNoZeroStep = 0; 
  }
 
  GetChordFinder()->SetDeltaChord(deltaChord);
  return TruePathLength;
}

// ---------------------------------------------------------------------------
// Dumps status of propagator
//
void
G4PropagatorInField::printStatus( const G4FieldTrack&      StartFT,
                                  const G4FieldTrack&      CurrentFT, 
                                        G4double           requestStep, 
                                        G4double           safety,
                                        G4int              stepNo, 
                                        G4VPhysicalVolume* startVolume)
{
  const G4int verboseLevel = fVerboseLevel;
  const G4ThreeVector StartPosition       = StartFT.GetPosition();
  const G4ThreeVector StartUnitVelocity   = StartFT.GetMomentumDir();
  const G4ThreeVector CurrentPosition     = CurrentFT.GetPosition();
  const G4ThreeVector CurrentUnitVelocity = CurrentFT.GetMomentumDir();

  G4double step_len = CurrentFT.GetCurveLength() - StartFT.GetCurveLength();

  G4long oldprec;   // cout/cerr precision settings
      
  if( ((stepNo == 0) && (verboseLevel <3)) || (verboseLevel >= 3) )
  {
    oldprec = G4cout.precision(4);
    G4cout << std::setw( 5) << "Step#" 
           << std::setw(10) << "  s  " << " "
           << std::setw(10) << "X(mm)" << " "
           << std::setw(10) << "Y(mm)" << " "  
           << std::setw(10) << "Z(mm)" << " "
           << std::setw( 7) << " N_x " << " "
           << std::setw( 7) << " N_y " << " "
           << std::setw( 7) << " N_z " << " " ;
    G4cout << std::setw( 7) << " Delta|N|" << " "
           << std::setw( 9) << "StepLen" << " "  
           << std::setw(12) << "StartSafety" << " "  
           << std::setw( 9) << "PhsStep" << " ";  
    if( startVolume != nullptr )
      { G4cout << std::setw(18) << "NextVolume" << " "; }
    G4cout.precision(oldprec);
    G4cout << G4endl;
  }
  if((stepNo == 0) && (verboseLevel <=3))
  {
    // Recurse to print the start values
    //
    printStatus( StartFT, StartFT, -1.0, safety, -1, startVolume);
  }
  if( verboseLevel <= 3 )
  {
    if( stepNo >= 0)
      { G4cout << std::setw( 4) << stepNo << " "; }
    else
      { G4cout << std::setw( 5) << "Start" ; }
    oldprec = G4cout.precision(8);
    G4cout << std::setw(10) << CurrentFT.GetCurveLength() << " "; 
    G4cout.precision(8);
    G4cout << std::setw(10) << CurrentPosition.x() << " "
           << std::setw(10) << CurrentPosition.y() << " "
           << std::setw(10) << CurrentPosition.z() << " ";
    G4cout.precision(4);
    G4cout << std::setw( 7) << CurrentUnitVelocity.x() << " "
           << std::setw( 7) << CurrentUnitVelocity.y() << " "
           << std::setw( 7) << CurrentUnitVelocity.z() << " ";
    G4cout.precision(3); 
    G4cout << std::setw( 7)
           << CurrentFT.GetMomentum().mag()-StartFT.GetMomentum().mag() << " "; 
    G4cout << std::setw( 9) << step_len << " "; 
    G4cout << std::setw(12) << safety << " ";
    if( requestStep != -1.0 ) 
      { G4cout << std::setw( 9) << requestStep << " "; }
    else
      { G4cout << std::setw( 9) << "Init/NotKnown" << " "; }
    if( startVolume != 0)
      { G4cout << std::setw(12) << startVolume->GetName() << " "; }
    G4cout.precision(oldprec);
    G4cout << G4endl;
  }
  else // if( verboseLevel > 3 )
  {
    //  Multi-line output
      
    G4cout << "Step taken was " << step_len  
           << " out of PhysicalStep = " <<  requestStep << G4endl;
    G4cout << "Final safety is: " << safety << G4endl;
    G4cout << "Chord length = " << (CurrentPosition-StartPosition).mag()
           << G4endl;
    G4cout << G4endl; 
  }
}

// ---------------------------------------------------------------------------
// Prints Step diagnostics
//
void 
G4PropagatorInField::PrintStepLengthDiagnostic(
                          G4double CurrentProposedStepLength,
                          G4double decreaseFactor,
                          G4double stepTrial,
                    const G4FieldTrack& )
{
  G4long iprec= G4cout.precision(8); 
  G4cout << " " << std::setw(12) << " PiF: NoZeroStep " 
         << " " << std::setw(20) << " CurrentProposed len " 
         << " " << std::setw(18) << " Full_curvelen_last" 
         << " " << std::setw(18) << " last proposed len " 
         << " " << std::setw(18) << " decrease factor   " 
         << " " << std::setw(15) << " step trial  " 
         << G4endl;

  G4cout << " " << std::setw(10) << fNoZeroStep << "  "
         << " " << std::setw(20) << CurrentProposedStepLength
         << " " << std::setw(18) << fFull_CurveLen_of_LastAttempt
         << " " << std::setw(18) << fLast_ProposedStepLength 
         << " " << std::setw(18) << decreaseFactor
         << " " << std::setw(15) << stepTrial
         << G4endl;
  G4cout.precision( iprec );
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

  if (fpTrajectoryFilter != nullptr)
  {
    return fpTrajectoryFilter->GimmeThePointsAndForgetThem();
  }
  else
  {
    return nullptr;
  }
}

// ---------------------------------------------------------------------------
//
void 
G4PropagatorInField::SetTrajectoryFilter(G4VCurvedTrajectoryFilter* filter)
{
  fpTrajectoryFilter = filter;
}

// ---------------------------------------------------------------------------
//
void G4PropagatorInField::ClearPropagatorState()
{
  // Goal: Clear all memory of previous steps,  cached information

  fParticleIsLooping = false;
  fNoZeroStep = 0;

  fSetFieldMgr = false;  // Has field-manager been set for the current step?
  fEpsilonStep= 1.0e-5;  // Relative accuracy of current Step
  
  End_PointAndTangent= G4FieldTrack( G4ThreeVector(0.,0.,0.),
                                     G4ThreeVector(0.,0.,0.),
                                     0.0,0.0,0.0,0.0,0.0); 
  fFull_CurveLen_of_LastAttempt = -1; 
  fLast_ProposedStepLength = -1;

  fPreviousSftOrigin= G4ThreeVector(0.,0.,0.);
  fPreviousSafety= 0.0;

  fNewTrack = true;
}

// ---------------------------------------------------------------------------
//
G4FieldManager* G4PropagatorInField::
FindAndSetFieldManager( G4VPhysicalVolume* pCurrentPhysicalVolume )
{
  G4FieldManager* currentFieldMgr;

  currentFieldMgr = fDetectorFieldMgr;
  if( pCurrentPhysicalVolume != nullptr )
  {
     G4FieldManager *pRegionFieldMgr = nullptr, *localFieldMgr = nullptr;
     G4LogicalVolume* pLogicalVol = pCurrentPhysicalVolume->GetLogicalVolume();

     if( pLogicalVol != nullptr )
     { 
        // Value for Region, if any, overrides
        //
        G4Region*  pRegion = pLogicalVol->GetRegion();
        if( pRegion != nullptr )
        { 
           pRegionFieldMgr = pRegion->GetFieldManager();
           if( pRegionFieldMgr != nullptr )
           {
              currentFieldMgr= pRegionFieldMgr;
           }
        }

        // 'Local' Value from logical volume, if any, overrides
        //
        localFieldMgr = pLogicalVol->GetFieldManager();
        if ( localFieldMgr != nullptr )
        {
           currentFieldMgr = localFieldMgr;
        }
     }
  }
  fCurrentFieldMgr = currentFieldMgr;

  // Flag that field manager has been set
  //
  fSetFieldMgr = true;

  return currentFieldMgr;
}

// ---------------------------------------------------------------------------
//
G4int G4PropagatorInField::SetVerboseLevel( G4int level )
{
  G4int oldval = fVerboseLevel;
  fVerboseLevel = level;

  // Forward the verbose level 'reduced' to ChordFinder,
  // MagIntegratorDriver ... ? 
  //
  auto integrDriver = GetChordFinder()->GetIntegrationDriver(); 
  integrDriver->SetVerboseLevel( fVerboseLevel - 2 );
  G4cout << "Set Driver verbosity to " << fVerboseLevel - 2 << G4endl;

  return oldval;
}

// ---------------------------------------------------------------------------
//
void G4PropagatorInField::ReportLoopingParticle( G4int count,
                                                 G4double StepTaken,
                                                 G4double StepRequested,
                                                 const char* methodName,
                                                 G4ThreeVector momentumVec,
                                                 G4VPhysicalVolume* pPhysVol )
{
   std::ostringstream message;
   G4double fraction = StepTaken / StepRequested;
   message << " Unfinished integration of track (likely looping particle)  "
           << " of momentum " << momentumVec << " ( magnitude = "
           << momentumVec.mag() << " ) " << G4endl
           << " after " << count << " field substeps "
           << " totaling " << std::setprecision(12) << StepTaken / mm << " mm "
           << " out of requested step " << std::setprecision(12)
           << StepRequested / mm << " mm ";
   message << " a fraction of ";
   G4int prec = 4;
   if( fraction > 0.99 )
   {
     prec = 7;
   }
   else
   {
     if (fraction > 0.97 )  { prec = 5; }
   }
   message << std::setprecision(prec) 
           << 100. * StepTaken / StepRequested << " % " << G4endl ;
   if( pPhysVol )
   {
      message << " in volume " << pPhysVol->GetName() ;
      auto material = pPhysVol->GetLogicalVolume()->GetMaterial();
      if( material != nullptr )
         message << " with material " << material->GetName()
                 << " ( density = "
                 << material->GetDensity() / ( g/(cm*cm*cm) ) << " g / cm^3 ) ";
   }
   else
   {
      message << " in unknown (null) volume. " ;
   }
   G4Exception(methodName, "GeomNav1002", JustWarning, message);   
}

// ---------------------------------------------------------------------------
//
void G4PropagatorInField::ReportStuckParticle( G4int    noZeroSteps,
                                               G4double proposedStep,
                                               G4double lastTriedStep,
                                               G4VPhysicalVolume* physVol )
{
   std::ostringstream message;
   message << "Particle is stuck; it will be killed." << G4endl
           << "  Zero progress for " << noZeroSteps << " attempted steps." 
           << G4endl
           << "  Proposed Step is " << proposedStep
           << " but Step Taken is "<< lastTriedStep << G4endl;
   if( physVol != nullptr )
      message << " in volume " << physVol->GetName() ; 
   else
      message << " in unknown or null volume. " ; 
   G4Exception("G4PropagatorInField::ComputeStep()",
               "GeomNav1002", JustWarning, message);
}
