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
// $Id: G4SimpleLocator.cc 96743 2016-05-03 08:01:33Z gcosmo $
//
// Class G4SimpleLocator implementation
//
// 27.10.08 - Tatiana Nikitina, extracted from G4PropagatorInField class
// 04.10.11 - John Apostolakis, revised convergence to use Surface Normal
// ---------------------------------------------------------------------------

#include <iomanip>

#include "G4ios.hh"
#include "G4SimpleLocator.hh"

G4SimpleLocator::G4SimpleLocator(G4Navigator *theNavigator)
  : G4VIntersectionLocator(theNavigator)
{
}      

G4SimpleLocator::~G4SimpleLocator()
{
}

// --------------------------------------------------------------------------
// G4bool G4PropagatorInField::LocateIntersectionPoint( 
//        const G4FieldTrack&       CurveStartPointVelocity,   // A
//        const G4FieldTrack&       CurveEndPointVelocity,     // B
//        const G4ThreeVector&      TrialPoint,                // E
//              G4FieldTrack&       IntersectedOrRecalculated  // Output
//              G4bool&             recalculated )             // Out
// --------------------------------------------------------------------------
//
// Function that returns the intersection of the true path with the surface
// of the current volume (either the external one or the inner one with one
// of the daughters:
//
//     A = Initial point
//     B = another point 
//
// Both A and B are assumed to be on the true path:
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
// NOTE: implementation taken from G4PropagatorInField
//
G4bool G4SimpleLocator::EstimateIntersectionPoint( 
         const  G4FieldTrack&       CurveStartPointVelocity,     // A
         const  G4FieldTrack&       CurveEndPointVelocity,       // B
         const  G4ThreeVector&      TrialPoint,                  // E
                G4FieldTrack&       IntersectedOrRecalculatedFT, // Output
                G4bool&             recalculatedEndPoint,        // Out
                G4double            &fPreviousSafety,            //In/Out
                G4ThreeVector       &fPreviousSftOrigin)         //In/Out
{
  // Find Intersection Point ( A, B, E )  of true path AB - start at E.

  G4bool found_approximate_intersection = false;
  G4bool there_is_no_intersection       = false;
  
  G4FieldTrack  CurrentA_PointVelocity = CurveStartPointVelocity; 
  G4FieldTrack  CurrentB_PointVelocity = CurveEndPointVelocity;
  G4ThreeVector CurrentE_Point = TrialPoint;
  G4bool        validNormalAtE = false;
  G4ThreeVector NormalAtEntry;

  G4FieldTrack  ApproxIntersecPointV(CurveEndPointVelocity); // FT-Def-Construct
  G4double      NewSafety = 0.0;
  G4bool last_AF_intersection = false;
  G4bool final_section = true;  // Shows whether current section is last
                                // (i.e. B=full end)
  recalculatedEndPoint = false; 

  G4bool restoredFullEndpoint = false;

  G4int substep_no = 0;

  // Limits for substep number
  //
  const G4int max_substeps  = 100000000;  // Test 120  (old value 100 )
  const G4int warn_substeps = 1000;       //      100  

  // Statistics for substeps
  //
  static G4ThreadLocal G4int max_no_seen= -1; 

  NormalAtEntry = GetSurfaceNormal( CurrentE_Point, validNormalAtE); 

#ifdef G4DEBUG_FIELD
  const G4double tolerance = 1.0e-8; 
  G4ThreeVector  StartPosition= CurveStartPointVelocity.GetPosition(); 
  if( (TrialPoint - StartPosition).mag() < tolerance * CLHEP::mm ) 
  {
     G4Exception("G4SimpleLocator::EstimateIntersectionPoint()", 
                 "GeomNav1002", JustWarning,
                 "Intersection point F is exactly at start point A." ); 
  }
#endif

  do 
  {
     G4ThreeVector Point_A = CurrentA_PointVelocity.GetPosition();  
     G4ThreeVector Point_B = CurrentB_PointVelocity.GetPosition();
       
     // F = a point on true AB path close to point E 
     // (the closest if possible)
     //
     ApproxIntersecPointV = GetChordFinderFor()
                            ->ApproxCurvePointV( CurrentA_PointVelocity, 
                                                 CurrentB_PointVelocity, 
                                                 CurrentE_Point,
                                                 GetEpsilonStepFor());
       //  The above method is the key & most intuitive part ...
      
#ifdef G4DEBUG_FIELD
     if( ApproxIntersecPointV.GetCurveLength() > 
         CurrentB_PointVelocity.GetCurveLength() * (1.0 + tolerance) )
     {
       G4Exception("G4SimpleLocator::EstimateIntersectionPoint()", 
                   "GeomNav0003", FatalException,
                   "Intermediate F point is past end B point!" ); 
     }
#endif

     G4ThreeVector CurrentF_Point= ApproxIntersecPointV.GetPosition();

     // First check whether EF is small - then F is a good approx. point 
     // Calculate the length and direction of the chord AF
     //
     G4ThreeVector  ChordEF_Vector = CurrentF_Point - CurrentE_Point;

     G4ThreeVector  NewMomentumDir= ApproxIntersecPointV.GetMomentumDir(); 
     G4double       MomDir_dot_Norm= NewMomentumDir.dot( NormalAtEntry ) ;

     G4ThreeVector  ChordAB           = Point_B - Point_A;

#ifdef G4DEBUG_FIELD
     G4VIntersectionLocator::
       ReportTrialStep( substep_no, ChordAB, ChordEF_Vector, 
                      NewMomentumDir, NormalAtEntry, validNormalAtE ); 
#endif
     // Check Sign is always exiting !! TODO
     // Could ( > -epsilon) be used instead?
     //
     G4bool adequate_angle = ( MomDir_dot_Norm >= 0.0 ) 
                          || (! validNormalAtE );  // Invalid
     G4double EF_dist2= ChordEF_Vector.mag2();
     if ( ( EF_dist2  <= sqr(fiDeltaIntersection) && ( adequate_angle ) )
       || ( EF_dist2 <= kCarTolerance*kCarTolerance ) )
     {
       found_approximate_intersection = true;
        
       // Create the "point" return value
       //
       IntersectedOrRecalculatedFT = ApproxIntersecPointV;
       IntersectedOrRecalculatedFT.SetPosition( CurrentE_Point );

       if ( GetAdjustementOfFoundIntersection() )
       {
         // Try to Get Correction of IntersectionPoint using SurfaceNormal()
         //  
         G4ThreeVector IP;
         G4ThreeVector MomentumDir= ApproxIntersecPointV.GetMomentumDirection();
         G4bool goodCorrection = AdjustmentOfFoundIntersection( Point_A,
                                   CurrentE_Point, CurrentF_Point, MomentumDir,
                                   last_AF_intersection, IP, NewSafety,
                                   fPreviousSafety, fPreviousSftOrigin );

         if(goodCorrection)
         {
           IntersectedOrRecalculatedFT = ApproxIntersecPointV;
           IntersectedOrRecalculatedFT.SetPosition(IP);
         }
       }

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
       GetNavigatorFor()->LocateGlobalPointWithinVolume( Point_A );

       G4ThreeVector PointG;   // Candidate intersection point
       G4double stepLengthAF; 
       G4bool usedNavigatorAF = false; 
       G4bool Intersects_AF = IntersectChord( Point_A,   
                                              CurrentF_Point,
                                              NewSafety,
                                              fPreviousSafety,
                                              fPreviousSftOrigin,
                                              stepLengthAF,
                                              PointG,
                                              &usedNavigatorAF );
       last_AF_intersection = Intersects_AF;
       if( Intersects_AF )
       {
         // G is our new Candidate for the intersection point.
         // It replaces  "E" and we will repeat the test to see if
         // it is a good enough approximate point for us.
         //       B    <- F
         //       E    <- G

         CurrentB_PointVelocity = ApproxIntersecPointV;
         CurrentE_Point = PointG;

         // Need to recalculate the Exit Normal at the new PointG 
         // Relies on a call to Navigator::ComputeStep in IntersectChord above
         // If the safety was adequate (for the step) this would NOT be called!
         // But this must not occur, no intersection can be found in that case,
         // so this branch, ie if( Intersects_AF ) would not be reached!
         //
         G4bool validNormalLast; 
         NormalAtEntry  = GetSurfaceNormal( PointG, validNormalLast ); 
         validNormalAtE = validNormalLast; 

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
         GetNavigatorFor()->LocateGlobalPointWithinVolume( CurrentF_Point );

         G4double stepLengthFB;
         G4ThreeVector PointH;
         G4bool usedNavigatorFB=false; 

         // Check whether any volumes are encountered by the chord FB
         // ---------------------------------------------------------

         G4bool Intersects_FB = IntersectChord( CurrentF_Point, Point_B, 
                                                NewSafety,fPreviousSafety,
                                                fPreviousSftOrigin,
                                                stepLengthFB,
                                                PointH, &usedNavigatorFB );
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
           //
           CurrentA_PointVelocity = ApproxIntersecPointV;
           CurrentE_Point = PointH;

           // Need to recalculate the Exit Normal at the new PointG
           // Relies on call to Navigator::ComputeStep in IntersectChord above
           // If safety was adequate (for the step) this would NOT be called!
           // But this must not occur, no intersection found in that case,
           // so this branch, i.e. if( Intersects_AF ) would not be reached!
           //
           G4bool validNormalLast; 
           NormalAtEntry  = GetSurfaceNormal( PointH, validNormalLast ); 
           validNormalAtE = validNormalLast;
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
              
             there_is_no_intersection = true;   // real final_section
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
       if( curveDist*curveDist*(1+2* GetEpsilonStepFor()) < linDistSq )
       {
         // Re-integrate to obtain a new B
         //
         G4FieldTrack newEndPointFT =
                 ReEstimateEndpoint( CurrentA_PointVelocity,
                                     CurrentB_PointVelocity,
                                     linDistSq,    // to avoid recalculation
                                     curveDist );
         G4FieldTrack oldPointVelB = CurrentB_PointVelocity; 
         CurrentB_PointVelocity = newEndPointFT;

         if( (final_section)) // real final section
         {
           recalculatedEndPoint = true;
           IntersectedOrRecalculatedFT = newEndPointFT;
             // So that we can return it, if it is the endpoint!
         }
       }
       if( curveDist < 0.0 )
       {
         fVerboseLevel = 5; // Print out a maximum of information
         printStatus( CurrentA_PointVelocity,  CurrentB_PointVelocity,
                      -1.0, NewSafety,  substep_no );
         std::ostringstream message;
         message << "Error in advancing propagation." << G4endl
                 << "        Point A (start) is " << CurrentA_PointVelocity
                 << G4endl
                 << "        Point B (end)   is " << CurrentB_PointVelocity
                 << G4endl
                 << "        Curve distance is " << curveDist << G4endl
                 << G4endl
                 << "The final curve point is not further along"
                 << " than the original!" << G4endl;

         if( recalculatedEndPoint )
         {
           message << "Recalculation of EndPoint was called with fEpsStep= "
                   << GetEpsilonStepFor() << G4endl;
         }
         message.precision(20);
         message << " Point A (Curve start)   is " << CurveStartPointVelocity
                 << G4endl
                 << " Point B (Curve   end)   is " << CurveEndPointVelocity
                 << G4endl
                 << " Point A (Current start) is " << CurrentA_PointVelocity
                 << G4endl
                 << " Point B (Current end)   is " << CurrentB_PointVelocity
                 << G4endl
                 << " Point E (Trial Point)   is " << CurrentE_Point
                 << G4endl
                 << " Point F (Intersection)  is " << ApproxIntersecPointV
                 << G4endl
                 << "        LocateIntersection parameters are : Substep no= "
                 << substep_no;

         G4Exception("G4SimpleLocator::EstimateIntersectionPoint()",
                     "GeomNav0003", FatalException, message);
       }

       if(restoredFullEndpoint)
       {
         final_section = restoredFullEndpoint;
         restoredFullEndpoint = false;
       }
     } // EndIf ( E is close enough to the curve, ie point F. )
       // tests ChordAF_Vector.mag() <= maximum_lateral_displacement 

#ifdef G4DEBUG_LOCATE_INTERSECTION  
     G4int trigger_substepno_print= warn_substeps - 20;

     if( substep_no >= trigger_substepno_print )
     {
       G4cout << "Difficulty in converging in "
              << "G4SimpleLocator::EstimateIntersectionPoint():"
              << G4endl
              << "    Substep no = " << substep_no << G4endl;
       if( substep_no == trigger_substepno_print )
       {
         printStatus( CurveStartPointVelocity, CurveEndPointVelocity,
                      -1.0, NewSafety, 0);
       }
       G4cout << " State of point A: "; 
       printStatus( CurrentA_PointVelocity, CurrentA_PointVelocity,
                    -1.0, NewSafety, substep_no-1, 0);
       G4cout << " State of point B: "; 
       printStatus( CurrentA_PointVelocity, CurrentB_PointVelocity,
                    -1.0, NewSafety, substep_no);
     }
#endif
     substep_no++; 

  } while ( ( ! found_approximate_intersection )
           && ( ! there_is_no_intersection )     
           && ( substep_no <= max_substeps) ); // UNTIL found or failed

  if( substep_no > max_no_seen )
  {
    max_no_seen = substep_no; 
#ifdef G4DEBUG_LOCATE_INTERSECTION  
    if( max_no_seen > warn_substeps )
    {
      trigger_substepno_print = max_no_seen-20; // Want to see last 20 steps 
    } 
#endif
  }

  if(  ( substep_no >= max_substeps)
      && !there_is_no_intersection
      && !found_approximate_intersection )
  {
    G4cout << "ERROR - G4SimpleLocator::EstimateIntersectionPoint()" << G4endl
           << "        Start and Endpoint of Requested Step:" << G4endl;
    printStatus( CurveStartPointVelocity, CurveEndPointVelocity,
                 -1.0, NewSafety, 0);
    G4cout << G4endl
           << "        Start and end-point of current Sub-Step:" << G4endl;
    printStatus( CurrentA_PointVelocity, CurrentA_PointVelocity,
                 -1.0, NewSafety, substep_no-1);
    printStatus( CurrentA_PointVelocity, CurrentB_PointVelocity,
                 -1.0, NewSafety, substep_no);

    std::ostringstream message;
    message << "Convergence is requiring too many substeps: "
            << substep_no << G4endl
            << "          Abandoning effort to intersect." << G4endl
            << "          Found intersection = "
            << found_approximate_intersection << G4endl
            << "          Intersection exists = "
            << !there_is_no_intersection << G4endl;
    message.precision(10); 
    G4double done_len = CurrentA_PointVelocity.GetCurveLength(); 
    G4double full_len = CurveEndPointVelocity.GetCurveLength();
    message << "          Undertaken only length: " << done_len
            << " out of " << full_len << " required." << G4endl
            << "          Remaining length = " << full_len-done_len; 

    G4Exception("G4SimpleLocator::EstimateIntersectionPoint()",
                "GeomNav0003", FatalException, message);
  }
  else if( substep_no >= warn_substeps )
  {  
    std::ostringstream message;
    message.precision(10); 

    message << "Many substeps while trying to locate intersection." << G4endl
            << "          Undertaken length: "  
            << CurrentB_PointVelocity.GetCurveLength() 
            << " - Needed: "  << substep_no << " substeps." << G4endl
            << "          Warning level = " << warn_substeps
            << " and maximum substeps = " << max_substeps;
    G4Exception("G4SimpleLocator::EstimateIntersectionPoint()",
                "GeomNav1002", JustWarning, message);
  }
  return  !there_is_no_intersection; //  Success or failure
}
