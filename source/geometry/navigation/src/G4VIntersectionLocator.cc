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
// $Id: G4VIntersectionLocator.cc,v 1.4.2.1 2009/08/11 13:31:19 gcosmo Exp $
// GEANT4 tag $Name: geant4-09-02-patch-02 $
//
// Class G4VIntersectionLocator implementation
//
// 27.10.08 - John Apostolakis, Tatiana Nikitina.
// ---------------------------------------------------------------------------
 
#include <iomanip>

#include "globals.hh"
#include "G4ios.hh"
#include "G4VIntersectionLocator.hh"
#include "G4GeometryTolerance.hh"

///////////////////////////////////////////////////////////////////////////
//
// Constructor
//
G4VIntersectionLocator:: G4VIntersectionLocator(G4Navigator *theNavigator)
{
  kCarTolerance = G4GeometryTolerance::GetInstance()->GetSurfaceTolerance();
  fiNavigator = theNavigator;
  fVerboseLevel = 0;
  fUseNormalCorrection = false;
  fHelpingNavigator = new G4Navigator();
}      

///////////////////////////////////////////////////////////////////////////
//
// Destructor.
//
G4VIntersectionLocator::~G4VIntersectionLocator()
{
  delete fHelpingNavigator;
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
      
  if( ((stepNo == 0) && (verboseLevel <3)) || (verboseLevel >= 3) )
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
  
    G4cout << G4endl;
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
     G4cout << std::setw( 7)
            << CurrentFT.GetMomentum().mag()- StartFT.GetMomentum().mag()
            << " "; 
     //   << std::setw( 7)
     //   << CurrentUnitVelocity.z() - InitialUnitVelocity.z() << " ";
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
       advanceLength=(EstimatedEndStateB.GetPosition()
                     -newEndPoint.GetPosition()).mag();
     }
     goodAdvance= 
       integrDriver->AccurateAdvance(newEndPoint, advanceLength,
                                     GetEpsilonStepFor());
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

///////////////////////////////////////////////////////////////////////////
//
// Method for finding SurfaceNormal of Intersecting Solid 
//
G4ThreeVector G4VIntersectionLocator::
GetLocalSurfaceNormal(const G4ThreeVector &CurrentE_Point, G4bool &validNormal)
{
  G4ThreeVector Normal(G4ThreeVector(0,0,0));
  G4VPhysicalVolume* located;

  validNormal = false;
  fHelpingNavigator->SetWorldVolume(GetNavigatorFor()->GetWorldVolume());
  located = fHelpingNavigator->LocateGlobalPointAndSetup( CurrentE_Point );
  G4TouchableHistoryHandle aTouchable = fHelpingNavigator
                                      ->CreateTouchableHistoryHandle();
  G4ThreeVector localPosition = aTouchable->GetHistory()
                ->GetTopTransform().TransformPoint(CurrentE_Point);

  if( located != 0)
  { 
    if (located->GetLogicalVolume()
        ->GetSolid()->Inside(localPosition)==kSurface)
    {
      Normal = located->GetLogicalVolume()
                      ->GetSolid()->SurfaceNormal(localPosition);
      validNormal = true;
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
  Normal=GetLocalSurfaceNormal(CurrentE_Point,validNormal);
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
