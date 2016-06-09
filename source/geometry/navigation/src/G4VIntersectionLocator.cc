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
// $Id: G4VIntersectionLocator.cc,v 1.8 2010-07-13 15:59:42 gcosmo Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// Class G4VIntersectionLocator implementation
//
// 27.10.08 - John Apostolakis, Tatiana Nikitina.
// ---------------------------------------------------------------------------
 
#include <iomanip>
#include <sstream>

#include "globals.hh"
#include "G4ios.hh"
#include "G4VIntersectionLocator.hh"
#include "G4GeometryTolerance.hh"

///////////////////////////////////////////////////////////////////////////
//
// Constructor
//
G4VIntersectionLocator:: G4VIntersectionLocator(G4Navigator *theNavigator): 
  fUseNormalCorrection(false), 
  fiNavigator( theNavigator ),
  fiChordFinder( 0 ),            // Not set - overridden at each step
  fiEpsilonStep( -1.0 ),         // Out of range - overridden at each step
  fiDeltaIntersection( -1.0 ),   // Out of range - overridden at each step
  fiUseSafety(false),            // Default - overridden at each step
  fpTouchable(0)           
{
  kCarTolerance = G4GeometryTolerance::GetInstance()->GetSurfaceTolerance();
  fVerboseLevel = 0;
  fHelpingNavigator = new G4Navigator();
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
// Dumps status of propagator.
//
void
G4VIntersectionLocator::printStatus( const G4FieldTrack&        StartFT,
                                     const G4FieldTrack&        CurrentFT, 
                                           G4double             requestStep, 
                                           G4double             safety,
                                           G4int                stepNo)
{
  const G4int verboseLevel= fVerboseLevel;
  const G4ThreeVector StartPosition       = StartFT.GetPosition();
  const G4ThreeVector StartUnitVelocity   = StartFT.GetMomentumDir();
  const G4ThreeVector CurrentPosition     = CurrentFT.GetPosition();
  const G4ThreeVector CurrentUnitVelocity = CurrentFT.GetMomentumDir();

  G4double step_len = CurrentFT.GetCurveLength() - StartFT.GetCurveLength();
  G4int oldprc;  // cout/cerr precision settings

  if( ((stepNo == 0) && (verboseLevel <3)) || (verboseLevel >= 3) )
  {
    oldprc = G4cout.precision(4);
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
    G4cout << std::setw( 7) << " Delta|N|" << " "
           << std::setw( 9) << "StepLen" << " "  
           << std::setw(12) << "StartSafety" << " "  
           << std::setw( 9) << "PhsStep" << " ";  
    G4cout << G4endl;
    G4cout.precision(oldprc);
  }
  if((stepNo == 0) && (verboseLevel <=3))
  {
    // Recurse to print the start values
    //
    printStatus( StartFT, StartFT, -1.0, safety, -1);
  }
  if( verboseLevel <= 3 )
  {
    if( stepNo >= 0)
    {
       G4cout << std::setw( 4) << stepNo << " ";
    }
    else
    {
       G4cout << std::setw( 5) << "Start" ;
    }
    oldprc = G4cout.precision(8);
    G4cout << std::setw(10) << CurrentFT.GetCurveLength() << " "; 
    G4cout << std::setw(10) << CurrentPosition.x() << " "
           << std::setw(10) << CurrentPosition.y() << " "
           << std::setw(10) << CurrentPosition.z() << " ";
    G4cout.precision(4);
    G4cout << std::setw( 7) << CurrentUnitVelocity.x() << " "
           << std::setw( 7) << CurrentUnitVelocity.y() << " "
           << std::setw( 7) << CurrentUnitVelocity.z() << " ";
    G4cout.precision(3); 
    G4cout << std::setw( 7)
           << CurrentFT.GetMomentum().mag()- StartFT.GetMomentum().mag()
           << " "; 
    G4cout << std::setw( 9) << step_len << " "; 
    G4cout << std::setw(12) << safety << " ";
    if( requestStep != -1.0 )
    {
      G4cout << std::setw( 9) << requestStep << " ";
    }
    else
    {
      G4cout << std::setw( 9) << "Init/NotKnown" << " "; 
    }
    G4cout << G4endl;
    G4cout.precision(oldprc);
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
// ReEstimateEndPoint.
//
G4FieldTrack G4VIntersectionLocator::
ReEstimateEndpoint( const G4FieldTrack &CurrentStateA,  
                    const G4FieldTrack &EstimatedEndStateB,
                          G4double      linearDistSq,
                          G4double      curveDist )
{  
  G4FieldTrack newEndPoint( CurrentStateA );
  G4MagInt_Driver* integrDriver= GetChordFinderFor()->GetIntegrationDriver(); 

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
       goodAdvance=true;
     }
     else{
     goodAdvance= 
       integrDriver->AccurateAdvance(newEndPoint, advanceLength,
                                     GetEpsilonStepFor());
     }
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
  //  below are some diagnostics only -- before the return!
  // 
  static const G4String MethodName("G4VIntersectionLocator::ReEstimateEndpoint");

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
  if (  (noInaccuracyWarnings < maxNoWarnings ) || (fVerboseLevel > 1) )
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

///////////////////////////////////////////////////////////////////////////
//
// Method for finding SurfaceNormal of Intersecting Solid 
//
G4ThreeVector G4VIntersectionLocator::
GetLocalSurfaceNormal(const G4ThreeVector &CurrentE_Point, G4bool &validNormal)
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

    if( (pLogical != 0) && ( (pSolid=pLogical->GetSolid()) !=0 )  )
    {
      // G4bool     goodPoint,    nearbyPoint;   
      // G4int   numGoodPoints,   numNearbyPoints;  // --> use for stats
      if ( ( pSolid->Inside(localPosition)==kSurface )
           || ( pSolid->DistanceToOut(localPosition) < 1000.0 * kCarTolerance )
         )
      {
        Normal = pSolid->SurfaceNormal(localPosition);
        validNormal = true;

#ifdef G4DEBUG_FIELD
        if( std::fabs(Normal.mag2() - 1.0 ) > perMille) 
        {
          G4cerr << "PROBLEM in G4VIntersectionLocator::GetLocalSurfaceNormal."
                 << G4endl;
          G4cerr << "  Normal is not unit - mag=" << Normal.mag() << G4endl; 
          G4cerr << "  at trial local point " << CurrentE_Point << G4endl;
          G4cerr << "  Solid is " << *pSolid << G4endl;
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
AdjustmentOfFoundIntersection( const G4ThreeVector &CurrentA_Point,
                               const G4ThreeVector &CurrentE_Point, 
                               const G4ThreeVector &CurrentF_Point,
                               const G4ThreeVector &MomentumDir,
                               const G4bool         IntersectAF,
                                     G4ThreeVector &IntersectionPoint,  // I/O
                                     G4double      &NewSafety,          // I/O 
                                     G4double      &fPreviousSafety,    // I/O
                                     G4ThreeVector &fPreviousSftOrigin )// I/O
{     
  G4double dist,lambda;
  G4ThreeVector Normal, NewPoint, Point_G;
  G4bool goodAdjust=false, Intersects_FP=false, validNormal=false;

  // Get SurfaceNormal of Intersecting Solid
  //
  Normal = GetGlobalSurfaceNormal(CurrentE_Point,validNormal);
  if(!validNormal) { return false; }

  // Intersection between Line and Plane
  //
  G4double n_d_m = Normal.dot(MomentumDir);
  if ( std::abs(n_d_m)>kCarTolerance )
  {
    if ( fVerboseLevel>1 )
    {
      G4cerr << "WARNING - "
             << "G4VIntersectionLocator::AdjustementOfFoundIntersection()"
             << G4endl
             << "        No intersection. Parallels lines!" << G4endl;
        return false;
    }
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

G4ThreeVector
G4VIntersectionLocator::GetSurfaceNormal(const G4ThreeVector &CurrentInt_Point,
                                               G4bool &validNormal) // const
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
    message << "G4VIntersectionLocator::GetSurfaceNormal -- identified problem."
            << G4endl;
    message << "PROBLEM: Normal is not unit - magnitude = " << NormalAtEntryLast.mag()
            << G4endl; 
    message << "   at trial intersection point " << CurrentInt_Point << G4endl;
    message << "   Obtained from Get *Last* Surface Normal." << G4endl; 
    G4Exception("G4VIntersectionLocator::GetGlobalSurfaceNormal()",
                "InvalidNormal", JustWarning, message);
  }
#endif

  if( validNormalLast ) 
  {
    NormalAtEntry=NormalAtEntryLast;  
    validNormal  = validNormalLast; 
  }
  return NormalAtEntry; 
}

G4ThreeVector G4VIntersectionLocator::
GetGlobalSurfaceNormal(const G4ThreeVector &CurrentE_Point,
                             G4bool &validNormal)
{
  G4ThreeVector     localNormal=
      GetLocalSurfaceNormal( CurrentE_Point, validNormal );
  G4AffineTransform localToGlobal=           //  Must use the same Navigator !!
      fHelpingNavigator->GetLocalToGlobalTransform();
  G4ThreeVector     globalNormal =
    localToGlobal.TransformAxis( localNormal );

#ifdef G4DEBUG_FIELD
  if( validNormal && ( std::fabs(globalNormal.mag2() - 1.0) > perThousand ) ) 
  {
    std::ostringstream message; 
    message << "****************************************************************"
            << G4endl;
    message << " Bad Normal in G4VIntersectionLocator::GetGlobalSurfaceNormal"
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
    message << "****************************************************************"
            << G4endl;
    G4Exception("G4VIntersectionLocator::GetGlobalSurfaceNormal()",
                "InvalidNormal", JustWarning, message);
  }
#endif

  return globalNormal;
}

G4ThreeVector 
G4VIntersectionLocator::GetLastSurfaceNormal( G4ThreeVector intersectPoint,
                                              G4bool &normalIsValid)   const
{
  G4ThreeVector normalVec;
  G4bool        validNorm;
  normalVec    = fiNavigator->GetGlobalExitNormal( intersectPoint, &validNorm ); 
  normalIsValid= validNorm;

  return normalVec;
}

void G4VIntersectionLocator::ReportTrialStep( G4int step_no, 
                                        const G4ThreeVector& ChordAB_v,
                                        const G4ThreeVector& ChordEF_v,
                                        const G4ThreeVector& NewMomentumDir,
                                        const G4ThreeVector& NormalAtEntry,
                                              G4bool validNormal )
{
  G4double       ABchord_length  = ChordAB_v.mag(); 
  G4double       MomDir_dot_Norm = NewMomentumDir.dot( NormalAtEntry ) ;
  G4double       MomDir_dot_ABchord;
  MomDir_dot_ABchord= (1.0 / ABchord_length) * NewMomentumDir.dot( ChordAB_v );

  std::ostringstream  outStream; 
  outStream // G4cout 
    << std::setw(6)  << " Step# "
    << std::setw(17) << " |ChordEF|(mag)" << "  "
    << std::setw(18) << " uMomentum.Normal" << "  "
    << std::setw(18) << " uMomentum.ABdir " << "  " 
    << std::setw(16) << " AB-dist         " << " " 
    << " Chord Vector (EF) " 
    << G4endl;
  outStream.precision(7); 
  outStream  // G4cout 
    << " " << std::setw(5)  << step_no           
    << " " << std::setw(18) << ChordEF_v.mag() 
    << " " << std::setw(18) << MomDir_dot_Norm    
    << " " << std::setw(18) << MomDir_dot_ABchord 
    << " " << std::setw(12) << ABchord_length     
    << " " << ChordEF_v
    << G4endl;
  outStream  // G4cout
    << " MomentumDir= " << " " << NewMomentumDir 
    << " Normal at Entry E= " << NormalAtEntry
    << " AB chord =   " << ChordAB_v
    << G4endl;
  G4cout << outStream.str();  // ostr_verbose;

  if( ( std::fabs(NormalAtEntry.mag2() - 1.0) > perThousand ) ) 
  {
    G4cerr << " PROBLEM in G4VIntersectionLocator::ReportTrialStep " << G4endl
           << "         Normal is not unit - mag=" <<  NormalAtEntry.mag() 
           << "         ValidNormalAtE = " << validNormal
           << G4endl; 
  }
  return; 
}
