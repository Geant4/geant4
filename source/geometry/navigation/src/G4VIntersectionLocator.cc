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
// Class G4VIntersectionLocator implementation
//
// 27.10.08 - John Apostolakis, Tatiana Nikitina.
// ---------------------------------------------------------------------------
 
#include <iomanip>
#include <sstream>

#include "globals.hh"
#include "G4ios.hh"
#include "G4AutoDelete.hh"
#include "G4SystemOfUnits.hh"
#include "G4VIntersectionLocator.hh"
#include "G4GeometryTolerance.hh"

///////////////////////////////////////////////////////////////////////////
//
// Constructor
//
G4VIntersectionLocator::G4VIntersectionLocator(G4Navigator* theNavigator)
 : fiNavigator(theNavigator)
{
  kCarTolerance = G4GeometryTolerance::GetInstance()->GetSurfaceTolerance();

  if( fiNavigator->GetExternalNavigation() == nullptr )
  {
    fHelpingNavigator = new G4Navigator();
  }
  else // Must clone the navigator, together with External Navigation
  {
    fHelpingNavigator = fiNavigator->Clone();
  }
}      

///////////////////////////////////////////////////////////////////////////
//
// Destructor.
//
G4VIntersectionLocator::~G4VIntersectionLocator()
{
  delete fHelpingNavigator;
  delete fpTouchable;
}

///////////////////////////////////////////////////////////////////////////
//
// Dump status of propagator to cout (old method)
//
void
G4VIntersectionLocator::printStatus( const G4FieldTrack& StartFT,
                                     const G4FieldTrack& CurrentFT, 
                                           G4double      requestStep, 
                                           G4double      safety,
                                           G4int         stepNo)
{
  std::ostringstream os; 
  printStatus( StartFT,CurrentFT,requestStep,safety,stepNo,os,fVerboseLevel);
  G4cout << os.str();
}

///////////////////////////////////////////////////////////////////////////
//
// Dumps status of propagator.
//
void
G4VIntersectionLocator::printStatus( const G4FieldTrack& StartFT,
                                     const G4FieldTrack& CurrentFT, 
                                           G4double      requestStep, 
                                           G4double      safety,
                                           G4int         stepNo,
                                           std::ostream& os,
                                           G4int         verboseLevel)
{
  // const G4int verboseLevel= fVerboseLevel;
  const G4ThreeVector StartPosition       = StartFT.GetPosition();
  const G4ThreeVector StartUnitVelocity   = StartFT.GetMomentumDir();
  const G4ThreeVector CurrentPosition     = CurrentFT.GetPosition();
  const G4ThreeVector CurrentUnitVelocity = CurrentFT.GetMomentumDir();

  G4double step_len = CurrentFT.GetCurveLength() - StartFT.GetCurveLength();
  G4long oldprc;  // cout/cerr precision settings

  if( ((stepNo == 0) && (verboseLevel <3)) || (verboseLevel >= 3) )
  {
    oldprc = os.precision(4);
    os << std::setw( 6)  << " " 
           << std::setw( 25) << " Current Position  and  Direction" << " "
           << G4endl; 
    os << std::setw( 5) << "Step#" 
           << std::setw(10) << "  s  " << " "
           << std::setw(10) << "X(mm)" << " "
           << std::setw(10) << "Y(mm)" << " "  
           << std::setw(10) << "Z(mm)" << " "
           << std::setw( 7) << " N_x " << " "
           << std::setw( 7) << " N_y " << " "
           << std::setw( 7) << " N_z " << " " ;
    os << std::setw( 7) << " Delta|N|" << " "
           << std::setw( 9) << "StepLen" << " "  
           << std::setw(12) << "StartSafety" << " "  
           << std::setw( 9) << "PhsStep" << " ";  
    os << G4endl;
    os.precision(oldprc);
  }
  if((stepNo == 0) && (verboseLevel <=3))
  {
    // Recurse to print the start values
    //
    printStatus( StartFT, StartFT, -1.0, safety, -1, os, verboseLevel);
  }
  if( verboseLevel <= 3 )
  {
    if( stepNo >= 0)
    {
       os << std::setw( 4) << stepNo << " ";
    }
    else
    {
       os << std::setw( 5) << "Start" ;
    }
    oldprc = os.precision(8);
    os << std::setw(10) << CurrentFT.GetCurveLength() << " "; 
    os << std::setw(10) << CurrentPosition.x() << " "
           << std::setw(10) << CurrentPosition.y() << " "
           << std::setw(10) << CurrentPosition.z() << " ";
    os.precision(4);
    os << std::setw( 7) << CurrentUnitVelocity.x() << " "
           << std::setw( 7) << CurrentUnitVelocity.y() << " "
           << std::setw( 7) << CurrentUnitVelocity.z() << " ";
    os.precision(3); 
    os << std::setw( 7)
           << CurrentFT.GetMomentum().mag()- StartFT.GetMomentum().mag()
           << " "; 
    os << std::setw( 9) << step_len << " "; 
    os << std::setw(12) << safety << " ";
    if( requestStep != -1.0 )
    {
      os << std::setw( 9) << requestStep << " ";
    }
    else
    {
      os << std::setw( 9) << "Init/NotKnown" << " "; 
    }
    os << G4endl;
    os.precision(oldprc);
  }
  else // if( verboseLevel > 3 )
  {
    //  Multi-line output
       
    os << "Step taken was " << step_len  
           << " out of PhysicalStep= " <<  requestStep << G4endl;
    os << "Final safety is: " << safety << G4endl;
    os << "Chord length = " << (CurrentPosition-StartPosition).mag()
           << G4endl;
    os << G4endl; 
  }
}

///////////////////////////////////////////////////////////////////////////
//
// ReEstimateEndPoint.
//
G4FieldTrack G4VIntersectionLocator::
ReEstimateEndpoint( const G4FieldTrack& CurrentStateA,  
                    const G4FieldTrack& EstimatedEndStateB,
                          G4double      , // linearDistSq,  // NOT used
                          G4double
#ifdef G4DEBUG_FIELD
  curveDist
#endif
                                   )
{  
  G4FieldTrack newEndPoint( CurrentStateA );
  auto integrDriver = GetChordFinderFor()->GetIntegrationDriver(); 

  G4FieldTrack retEndPoint( CurrentStateA );
  G4bool goodAdvance;
  G4int  itrial = 0;
  const G4int no_trials = 20;


  G4double endCurveLen= EstimatedEndStateB.GetCurveLength();

  do  // Loop checking, 07.10.2016, JA
  {
    G4double currentCurveLen = newEndPoint.GetCurveLength();
    G4double advanceLength = endCurveLen - currentCurveLen ; 
    if (std::abs(advanceLength)<kCarTolerance)
    {
      goodAdvance=true;
    }
    else
    {
      goodAdvance = integrDriver->AccurateAdvance(newEndPoint, advanceLength,
                                                  GetEpsilonStepFor());
    }
  }
  while( !goodAdvance && (++itrial < no_trials) );

  if( goodAdvance )
  {
    retEndPoint = newEndPoint; 
  }
  else
  {
    retEndPoint = EstimatedEndStateB; // Could not improve without major work !!
  }

  //  All the work is done
  //  below are some diagnostics only -- before the return!
  // 
  const G4String MethodName("G4VIntersectionLocator::ReEstimateEndpoint()");

#ifdef G4VERBOSE
  G4int  latest_good_trials = 0;
  if( itrial > 1)
  {
    if( fVerboseLevel > 0 )
    {
      G4cout << MethodName << " called - goodAdv= " << goodAdvance
             << " trials = " << itrial
             << " previous good= " << latest_good_trials
             << G4endl;
    }
    latest_good_trials = 0; 
  }
  else
  {   
    ++latest_good_trials; 
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
  G4double linearDist = ( EstimatedEndStateB.GetPosition() 
                      - CurrentStateA.GetPosition() ).mag(); 
  static G4int noInaccuracyWarnings = 0; 
  G4int maxNoWarnings = 10;
  if (  (noInaccuracyWarnings < maxNoWarnings ) 
       || (fVerboseLevel > 1) )
    {
      G4ThreeVector move = newEndPoint.GetPosition()
                         - EstimatedEndStateB.GetPosition();
      std::ostringstream message;
      message.precision(12);
      message << " Integration inaccuracy requires" 
              << " an adjustment in the step's endpoint."  << G4endl
              << "   Two mid-points are further apart than their"
              << " curve length difference"                << G4endl 
              << "   Dist = "       << linearDist
              << " curve length = " << curveDist             << G4endl; 
      message << " Correction applied is " << move.mag() << G4endl
              << "  Old Estimated B position= "
              << EstimatedEndStateB.GetPosition() << G4endl
              << "  Recalculated    Position= "
              << newEndPoint.GetPosition() << G4endl
              << "   Change ( new - old )   = " << move;
      G4Exception("G4VIntersectionLocator::ReEstimateEndpoint()",
                  "GeomNav1002", JustWarning, message);
    }
/*
#else
  // Statistics on the RMS value of the corrections

  static G4ThreadLocal G4int noCorrections = 0;
  ++noCorrections; 
  if( goodAdvance )
  {
    static G4ThreadLocal G4double sumCorrectionsSq;
    sumCorrectionsSq += (EstimatedEndStateB.GetPosition() - 
                         newEndPoint.GetPosition()).mag2();
  }
*/
#endif

  return retEndPoint;
}


///////////////////////////////////////////////////////////////////////////
//
// ReEstimateEndPoint.
//
//   The following values are returned in curveError 
//       0: Normal - no problem
//       1: Unexpected co-incidence - milder mixup
//       2: Real mixup - EndB is NOT past StartA  
//            ( ie. StartA's curve-lengh is bigger than EndB's)


G4bool G4VIntersectionLocator::
CheckAndReEstimateEndpoint( const G4FieldTrack& CurrentStartA,  
                            const G4FieldTrack& EstimatedEndB,
                                  G4FieldTrack& RevisedEndPoint,
                                  G4int&        curveError)
{
  G4double linDistSq, curveDist; 

  G4bool recalculated = false;
  curveError= 0;

  linDistSq = ( EstimatedEndB.GetPosition() 
              - CurrentStartA.GetPosition() ).mag2(); 
  curveDist = EstimatedEndB.GetCurveLength()
            - CurrentStartA.GetCurveLength();
  if(  (curveDist>=0.0) 
     && (curveDist*curveDist *(1.0+2.0*fiEpsilonStep ) < linDistSq ) )
  {
    G4FieldTrack newEndPointFT = EstimatedEndB;  // Unused

    if (curveDist>0.0) 
    {
       // Re-integrate to obtain a new B
       RevisedEndPoint = ReEstimateEndpoint( CurrentStartA,
                                             EstimatedEndB,
                                             linDistSq,    
                                             curveDist );
       recalculated = true;
    }
    else
    {
       // Zero length -> no advance!
       newEndPointFT = CurrentStartA;
       recalculated = true;
       curveError = 1;  // Unexpected co-incidence - milder mixup

       G4Exception("G4MultiLevelLocator::EstimateIntersectionPoint()",
           "GeomNav1002", JustWarning, 
           "A & B are at equal distance in 2nd half. A & B will coincide." );
    }
  }

  // Sanity check
  //
  if( curveDist < 0.0 )
  {
    curveError = 2;  //  Real mixup
  }
  return recalculated;
}

///////////////////////////////////////////////////////////////////////////
//
// Method for finding SurfaceNormal of Intersecting Solid 
//
G4ThreeVector G4VIntersectionLocator::
GetLocalSurfaceNormal(const G4ThreeVector& CurrentE_Point, G4bool& validNormal)
{
  G4ThreeVector Normal(G4ThreeVector(0.0,0.0,0.0));
  G4VPhysicalVolume* located;

  validNormal = false;
  fHelpingNavigator->SetWorldVolume(GetNavigatorFor()->GetWorldVolume());
  located = fHelpingNavigator->LocateGlobalPointAndSetup( CurrentE_Point );

  delete fpTouchable;
  fpTouchable = fHelpingNavigator->CreateTouchableHistory();

  // To check if we can use GetGlobalExitNormal() 
  //
  G4ThreeVector localPosition = fpTouchable->GetHistory()
                ->GetTopTransform().TransformPoint(CurrentE_Point);

  // Issue: in the case of coincident surfaces, this version does not recognise 
  //        which side you are located onto (can return vector with wrong sign.)
  // TO-DO: use direction (of chord) to identify volume we will be "entering"

  if( located != 0)
  { 
    G4LogicalVolume* pLogical= located->GetLogicalVolume(); 
    G4VSolid*        pSolid; 

    if( (pLogical != nullptr) && ( (pSolid=pLogical->GetSolid()) != nullptr ) )
    {
      if ( ( pSolid->Inside(localPosition)==kSurface )
        || ( pSolid->DistanceToOut(localPosition) < 1000.0 * kCarTolerance ) )
      {
        Normal = pSolid->SurfaceNormal(localPosition);
        validNormal = true;

#ifdef G4DEBUG_FIELD
        if( std::fabs(Normal.mag2() - 1.0 ) > CLHEP::perThousand) 
        {
          G4cerr << "PROBLEM in G4VIntersectionLocator::GetLocalSurfaceNormal."
                 << G4endl;
          G4cerr << "  Normal is not unit - mag=" << Normal.mag() << G4endl; 
          G4cerr << "  at trial local point " << CurrentE_Point << G4endl;
          G4cerr <<  "  Solid is " << *pSolid << G4endl;
        }
#endif
      }
    }
  }
  return Normal;
}

///////////////////////////////////////////////////////////////////////////
//
// Adjustment of Found Intersection
//
G4bool G4VIntersectionLocator::
AdjustmentOfFoundIntersection( const G4ThreeVector& CurrentA_Point,
                               const G4ThreeVector& CurrentE_Point, 
                               const G4ThreeVector& CurrentF_Point,
                               const G4ThreeVector& MomentumDir,
                               const G4bool         IntersectAF,
                                     G4ThreeVector& IntersectionPoint,  // I/O
                                     G4double&      NewSafety,          // I/O 
                                     G4double&      fPreviousSafety,    // I/O
                                     G4ThreeVector& fPreviousSftOrigin )// I/O
{     
  G4double dist,lambda;
  G4ThreeVector Normal, NewPoint, Point_G;
  G4bool goodAdjust = false, Intersects_FP = false, validNormal = false;

  // Get SurfaceNormal of Intersecting Solid
  //
  Normal = GetGlobalSurfaceNormal(CurrentE_Point,validNormal);
  if(!validNormal) { return false; }

  // Intersection between Line and Plane
  //
  G4double n_d_m = Normal.dot(MomentumDir);
  if ( std::abs(n_d_m)>kCarTolerance )
  {
#ifdef G4VERBOSE
    if ( fVerboseLevel>1 )
    {
      G4Exception("G4VIntersectionLocator::AdjustmentOfFoundIntersection()",
                  "GeomNav0003", JustWarning,
                  "No intersection. Parallels lines!");
    }
#endif
    lambda =- Normal.dot(CurrentF_Point-CurrentE_Point)/n_d_m;

    // New candidate for Intersection
    //
    NewPoint = CurrentF_Point+lambda*MomentumDir;

    // Distance from CurrentF to Calculated Intersection
    //
    dist = std::abs(lambda);

    if ( dist<kCarTolerance*0.001 )  { return false; }

    // Calculation of new intersection point on the path.
    //
    if ( IntersectAF )  //  First part intersects
    {
      G4double stepLengthFP; 
      G4ThreeVector Point_P = CurrentA_Point;
      GetNavigatorFor()->LocateGlobalPointWithinVolume(Point_P);
      Intersects_FP = IntersectChord( Point_P, NewPoint, NewSafety,
                                      fPreviousSafety, fPreviousSftOrigin,
                                      stepLengthFP, Point_G );

    }
    else   // Second part intersects
    {      
      G4double stepLengthFP; 
      GetNavigatorFor()->LocateGlobalPointWithinVolume(CurrentF_Point );
      Intersects_FP = IntersectChord( CurrentF_Point, NewPoint, NewSafety,
                                      fPreviousSafety, fPreviousSftOrigin,
                                      stepLengthFP, Point_G );
    }
    if ( Intersects_FP )
    {
      goodAdjust = true;
      IntersectionPoint = Point_G;              
    }
  }

  return goodAdjust;
}

///////////////////////////////////////////////////////////////////////////
//
// GetSurfaceNormal.
//
G4ThreeVector G4VIntersectionLocator::
GetSurfaceNormal(const G4ThreeVector& CurrentInt_Point,
                       G4bool& validNormal)
{
  G4ThreeVector NormalAtEntry; // ( -10. , -10., -10. ); 

  G4ThreeVector NormalAtEntryLast, NormalAtEntryGlobal, diffNormals;
  G4bool validNormalLast; 

  // Relies on a call to Navigator::ComputeStep in IntersectChord before
  // this call
  //
  NormalAtEntryLast = GetLastSurfaceNormal( CurrentInt_Point, validNormalLast );
    // May return valid=false in cases, including
    //  - if the candidate volume was not found (eg exiting world), or
    //  - a replica was involved -- determined the step size.
    // (This list is not complete.) 

#ifdef G4DEBUG_FIELD
  if  ( validNormalLast
   && ( std::fabs(NormalAtEntryLast.mag2() - 1.0) > perThousand ) )
  {
    std::ostringstream message; 
    message << "PROBLEM: Normal is not unit - magnitude = "
            << NormalAtEntryLast.mag() << G4endl; 
    message << "   at trial intersection point " << CurrentInt_Point << G4endl;
    message << "   Obtained from Get *Last* Surface Normal."; 
    G4Exception("G4VIntersectionLocator::GetSurfaceNormal()",
                "GeomNav1002", JustWarning, message);
  }
#endif

  if( validNormalLast ) 
  {
    NormalAtEntry = NormalAtEntryLast;  
  }
  validNormal = validNormalLast;

  return NormalAtEntry;
}

///////////////////////////////////////////////////////////////////////////
//
// GetGlobalSurfaceNormal.
//
G4ThreeVector G4VIntersectionLocator::
GetGlobalSurfaceNormal(const G4ThreeVector& CurrentE_Point,
                             G4bool& validNormal)
{
  G4ThreeVector localNormal = GetLocalSurfaceNormal(CurrentE_Point,validNormal);
  G4AffineTransform localToGlobal =          //  Must use the same Navigator !!
      fHelpingNavigator->GetLocalToGlobalTransform();
  G4ThreeVector globalNormal = localToGlobal.TransformAxis( localNormal );

#ifdef G4DEBUG_FIELD
  if( validNormal && ( std::fabs(globalNormal.mag2() - 1.0) > perThousand ) ) 
  {
    std::ostringstream message; 
    message << "**************************************************************"
            << G4endl;
    message << " Bad Normal in G4VIntersectionLocator::GetGlobalSurfaceNormal "
            << G4endl;
    message << "  * Constituents: " << G4endl;
    message << "    Local  Normal= " << localNormal << G4endl;
    message << "    Transform: " << G4endl
            << "      Net Translation= " << localToGlobal.NetTranslation()
            << G4endl
            << "      Net Rotation   = " << localToGlobal.NetRotation()
            << G4endl;
    message << "  * Result: " << G4endl;
    message << "     Global Normal= " << localNormal << G4endl;
    message << "**************************************************************";
    G4Exception("G4VIntersectionLocator::GetGlobalSurfaceNormal()",
                "GeomNav1002", JustWarning, message);
  }
#endif

  return globalNormal;
}

///////////////////////////////////////////////////////////////////////////
//
// GetLastSurfaceNormal.
//
G4ThreeVector G4VIntersectionLocator::
GetLastSurfaceNormal( const G4ThreeVector& intersectPoint,
                            G4bool& normalIsValid) const
{
  G4ThreeVector normalVec;
  G4bool validNorm;
  normalVec = fiNavigator->GetGlobalExitNormal( intersectPoint, &validNorm ); 
  normalIsValid = validNorm;

  return normalVec;
}

///////////////////////////////////////////////////////////////////////////
//
// ReportTrialStep.
//
void G4VIntersectionLocator::ReportTrialStep( G4int step_no, 
                                        const G4ThreeVector& ChordAB_v,
                                        const G4ThreeVector& ChordEF_v,
                                        const G4ThreeVector& NewMomentumDir,
                                        const G4ThreeVector& NormalAtEntry,
                                              G4bool validNormal )
{
  G4double ABchord_length  = ChordAB_v.mag(); 
  G4double MomDir_dot_Norm = NewMomentumDir.dot( NormalAtEntry );
  G4double MomDir_dot_ABchord;
  MomDir_dot_ABchord = (1.0 / ABchord_length) * NewMomentumDir.dot( ChordAB_v );

  std::ostringstream  outStream; 
  outStream << std::setw(6)  << " Step# "
    << std::setw(17) << " |ChordEF|(mag)" << "  "
    << std::setw(18) << " uMomentum.Normal" << "  "
    << std::setw(18) << " uMomentum.ABdir " << "  " 
    << std::setw(16) << " AB-dist         " << " " 
    << " Chord Vector (EF) " 
    << G4endl;
  outStream.precision(7); 
  outStream << " " << std::setw(5) << step_no           
    << " " << std::setw(18) << ChordEF_v.mag() 
    << " " << std::setw(18) << MomDir_dot_Norm    
    << " " << std::setw(18) << MomDir_dot_ABchord 
    << " " << std::setw(12) << ABchord_length     
    << " " << ChordEF_v
    << G4endl;
  outStream << " MomentumDir= " << " " << NewMomentumDir 
    << " Normal at Entry E= " << NormalAtEntry
    << " AB chord =   " << ChordAB_v
    << G4endl;
  G4cout << outStream.str();

  if( ( std::fabs(NormalAtEntry.mag2() - 1.0) > perThousand ) ) 
  {
    std::ostringstream message; 
    message << "Normal is not unit - mag= " << NormalAtEntry.mag() << G4endl
            << "         ValidNormalAtE = " << validNormal;
    G4Exception("G4VIntersectionLocator::ReportTrialStep()",
                "GeomNav1002", JustWarning, message);
  }
  return; 
}

///////////////////////////////////////////////////////////////////////////
//
// LocateGlobalPointWithinVolumeAndCheck.
//
// Locate point using navigator: updates state of Navigator
// By default, it assumes that the point is inside the current volume,
// and returns true.
// In check mode, checks whether the point is *inside* the volume.
//   If it is inside, it returns true
//   If not, issues a warning and returns false.
//
G4bool G4VIntersectionLocator::
LocateGlobalPointWithinVolumeAndCheck( const G4ThreeVector& position )
{
  G4bool good = true;
  G4Navigator* nav = GetNavigatorFor();
  const G4String
  MethodName("G4VIntersectionLocator::LocateGlobalPointWithinVolumeAndCheck()");

  if( fCheckMode )
  {
    G4bool navCheck= nav->IsCheckModeActive();  // Recover original value
    nav->CheckMode(true);

    // Identify the current volume
    
    G4TouchableHistoryHandle startTH= nav->CreateTouchableHistoryHandle();
    G4VPhysicalVolume* motherPhys = startTH->GetVolume();
    G4VSolid*          motherSolid = startTH->GetSolid();
    G4AffineTransform transform = nav->GetGlobalToLocalTransform();
    G4int motherCopyNo = motherPhys->GetCopyNo();
    
    // Let's check that the point is inside the current solid
    G4ThreeVector  localPosition = transform.TransformPoint(position);
    EInside        inMother = motherSolid->Inside( localPosition );
    if( inMother != kInside )
    {
      std::ostringstream message; 
      message << "Position located "
              << ( inMother == kSurface ? " on Surface " : " outside " )
              << "expected volume" << G4endl
              << "  Safety (from Outside) = "
              << motherSolid->DistanceToIn(localPosition);
      G4Exception("G4VIntersectionLocator::LocateGlobalPointWithinVolumeAndCheck()",
                  "GeomNav1002", JustWarning, message);
    }
    
    // 1. Simple next step - quick relocation and check result.
    // nav->LocateGlobalPointWithinVolume( position );
    
    // 2. Full relocation - to cross-check answer !
    G4VPhysicalVolume* nextPhysical= nav->LocateGlobalPointAndSetup(position);
    if(    (nextPhysical != motherPhys)
        || (nextPhysical->GetCopyNo() != motherCopyNo )
       )
    {
      G4Exception("G4VIntersectionLocator::LocateGlobalPointWithinVolumeAndCheck()",
                  "GeomNav1002", JustWarning,
                  "Position located outside expected volume.");
    }
    nav->CheckMode(navCheck);  // Recover original value
  }
  else
  {
    nav->LocateGlobalPointWithinVolume( position );
  }
  return good;
}

///////////////////////////////////////////////////////////////////////////
//
// LocateGlobalPointWithinVolumeCheckAndReport.
//
void G4VIntersectionLocator::
LocateGlobalPointWithinVolumeCheckAndReport( const G4ThreeVector& position,
                                             const G4String& CodeLocationInfo,
                                             G4int /* CheckMode */)
{
  // Save value of Check mode first
  G4bool oldCheck = GetCheckMode();
  
  G4bool ok = LocateGlobalPointWithinVolumeAndCheck( position );
  if( !ok )
  {
    std::ostringstream message; 
    message << "Failed point location." << G4endl
            << "   Code Location info: " << CodeLocationInfo;
    G4Exception("G4VIntersectionLocator::LocateGlobalPointWithinVolumeCheckAndReport()",
                "GeomNav1002", JustWarning, message);
  }
  
  SetCheckMode( oldCheck );
}

///////////////////////////////////////////////////////////////////////////
//
// ReportReversedPoints.
//
void G4VIntersectionLocator::
ReportReversedPoints( std::ostringstream& msg,
                      const G4FieldTrack& StartPointVel, 
                      const G4FieldTrack& EndPointVel,
                            G4double NewSafety, G4double epsStep,
                      const G4FieldTrack& A_PtVel,
                      const G4FieldTrack& B_PtVel,
                      const G4FieldTrack& SubStart_PtVel,
                      const G4ThreeVector& E_Point,
                      const G4FieldTrack& ApproxIntersecPointV,
                            G4int  substep_no, G4int substep_no_p, G4int depth )
{
   // Expect that 'msg' can hold the name of the calling method

   // FieldTrack 'points' A and B have been tangled
   // Whereas A should be before B, it is found that curveLen(B) < curveLen(A)
   G4int verboseLevel= 5; 
   G4double curveDist = B_PtVel.GetCurveLength() - A_PtVel.GetCurveLength();
   G4VIntersectionLocator::printStatus( A_PtVel,  B_PtVel,
                           -1.0, NewSafety,  substep_no, msg, verboseLevel );
   msg << "Error in advancing propagation." << G4endl
       << "   The final curve point is NOT further along"
       << "  than the original!" << G4endl
       << "   Going *backwards* from len(A) = " << A_PtVel.GetCurveLength()
       << "  to len(B) = " << B_PtVel.GetCurveLength() << G4endl
       << "      Curve distance is " << curveDist / CLHEP::millimeter << " mm "
       << G4endl
       << "      Point A' (start) is " << A_PtVel  << G4endl
       << "      Point B' (end)   is " << B_PtVel << G4endl;
   msg << "      fEpsStep= " << epsStep << G4endl << G4endl;

   G4long oldprc = msg.precision(20);
   msg << " In full precision, the position, momentum, E_kin, length, rest mass "
       << " ... are: " << G4endl;
   msg << " Point A[0] (Curve   start) is " << StartPointVel << G4endl
       << " Point S    (Sub     start) is " << SubStart_PtVel
       << " Point A'   (Current start) is " << A_PtVel << G4endl
       << " Point E    (Trial Point)   is " << E_Point << G4endl
       << " Point F    (Intersection)  is " << ApproxIntersecPointV << G4endl
       << " Point B'   (Current end)   is " << B_PtVel << G4endl
       << " Point B[0] (Curve   end)   is " << EndPointVel << G4endl
       << G4endl
       << " LocateIntersection parameters are : " << G4endl
       << "      Substep no (total) = "  << substep_no << G4endl
       << "      Substep no         = "  << substep_no_p << " at depth= " << depth;
   msg.precision(oldprc);
}

///////////////////////////////////////////////////////////////////////////
//
// ReportProgress.
//
void G4VIntersectionLocator::ReportProgress( std::ostream& oss,
                                    const G4FieldTrack& StartPointVel, 
                                    const G4FieldTrack& EndPointVel,
                                          G4int         substep_no, 
                                    const G4FieldTrack& A_PtVel,
                                    const G4FieldTrack& B_PtVel,
                                          G4double      safetyLast,
                                          G4int         depth )

{
  oss << "ReportProgress: Current status of intersection search: " << G4endl;
  if( depth > 0 ) oss << " Depth= " << depth;
  oss << " Substep no = " << substep_no << G4endl;
  G4int  verboseLevel = 5; 
  G4double safetyPrev = -1.0;  // Add as argument ?

  printStatus( StartPointVel, EndPointVel, -1.0, -1.0, -1, 
              oss, verboseLevel);  
  oss << " * Start and end-point of requested Step:" << G4endl;
  oss << " ** State of point A: "; 
  printStatus( A_PtVel, A_PtVel, -1.0, safetyPrev, substep_no-1,
               oss, verboseLevel);  
  oss << " ** State of point B: "; 
  printStatus( A_PtVel, B_PtVel, -1.0, safetyLast, substep_no, 
               oss, verboseLevel);
}

///////////////////////////////////////////////////////////////////////////
//
// ReportImmediateHit.
//
void
G4VIntersectionLocator::ReportImmediateHit( const char*          MethodName, 
                                            const G4ThreeVector& StartPosition, 
                                            const G4ThreeVector& TrialPoint, 
                                                  G4double       tolerance,
                                            unsigned long int    numCalls )
{
  static G4ThreadLocal unsigned int occurredOnTop= 0;
  static G4ThreadLocal G4ThreeVector* ptrLast = nullptr;
  if( ptrLast == nullptr )
  {
     ptrLast= new G4ThreeVector( DBL_MAX, DBL_MAX, DBL_MAX );
     G4AutoDelete::Register(ptrLast);
  }
  G4ThreeVector &lastStart= *ptrLast;

  if( (TrialPoint - StartPosition).mag2() < tolerance*tolerance) 
  {
     static G4ThreadLocal unsigned int numUnmoved = 0;
     static G4ThreadLocal unsigned int numStill = 0;    // Still at same point

     G4cout << "Intersection F == start A in " << MethodName;
     G4cout << "Start Point: " << StartPosition << G4endl;
     G4cout << " Start-Trial: " << TrialPoint - StartPosition;
     G4cout << " Start-last: " << StartPosition - lastStart;

     if( (StartPosition - lastStart).mag() < tolerance )
     {
        // We are at position of last 'Start' position - ie unmoved
        ++numUnmoved;
        ++numStill; 
        G4cout << " { Unmoved: "  << " still#= " << numStill
               << " total # = " << numUnmoved << " } - ";
     }
     else
     {
        numStill = 0; 
     }
     G4cout << " Occurred: " << ++occurredOnTop;  
     G4cout <<  " out of total calls= " << numCalls;
     G4cout << G4endl;
     lastStart = StartPosition;
  }
}  // End of ReportImmediateHit() 
