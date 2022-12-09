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
// Class G4MultiLevelLocator implementation
//
// 27.10.08 - Tatiana Nikitina.
// 04.10.11 - John Apostolakis, revised convergence to use Surface Normal
// ---------------------------------------------------------------------------

#include <iomanip>

#include "G4ios.hh"
#include "G4MultiLevelLocator.hh"
#include "G4LocatorChangeRecord.hh"
#include "G4LocatorChangeLogger.hh"

G4MultiLevelLocator::G4MultiLevelLocator(G4Navigator *theNavigator)
  : G4VIntersectionLocator(theNavigator)
{
  // In case of too slow progress in finding Intersection Point
  // intermediates Points on the Track must be stored.
  // Initialise the array of Pointers [max_depth+1] to do this  
  
  G4ThreeVector zeroV(0.0,0.0,0.0);
  for ( auto idepth=0; idepth<max_depth+1; ++idepth )
  {
    ptrInterMedFT[ idepth ] = new G4FieldTrack( zeroV, zeroV, 0., 0., 0., 0.);
  }

  if (fCheckMode)
  {
    //  Trial values       Loose         Medium      Tight
    //  To happen:         Infrequent    Occasional  Often
    SetMaxSteps(150);   //  300             85          25
    SetWarnSteps(80);   //  250             60          15
  }  
}      

G4MultiLevelLocator::~G4MultiLevelLocator()
{
  for ( auto idepth=0; idepth<max_depth+1; ++idepth )
  {
    delete ptrInterMedFT[idepth];
  }
#ifdef G4DEBUG_FIELD
  ReportStatistics(); 
#endif  
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
//        Potential reasons:
//        a) no segment found an intersection
//        b) too many steps were required - after that it abandoned the effort
//           and is returning how far it could go.  (New - 29 Oct 2015)
//           (If so, it must set 'recalculated' to true.) 
//        TODO/idea: add a new flag: 'unfinished' to identify these cases.
//
//        IntersectedOrRecalculated means different things:
//        a) if it is the same curve lenght along, it is a revision of the
//            original enpdoint due to the need for re-integration.
//        b) if it is at a shorter curve length, it is 'end of what it could do'
//           i.e. as far as it could go, because it took too many steps!
//        Note: IntersectedOrRecalculated is valid only if 'recalculated' is
//       'true'.
// --------------------------------------------------------------------------
// NOTE: implementation taken from G4PropagatorInField
//
G4bool G4MultiLevelLocator::EstimateIntersectionPoint( 
         const  G4FieldTrack&       CurveStartPointVelocity,       // A
         const  G4FieldTrack&       CurveEndPointVelocity,         // B
         const  G4ThreeVector&      TrialPoint,                    // E
                G4FieldTrack&       IntersectedOrRecalculatedFT,   // Output
                G4bool&             recalculatedEndPoint,          // Out
                G4double&           previousSafety,                // In/Out
                G4ThreeVector&      previousSftOrigin)             // In/Out
{
  // Find Intersection Point ( A, B, E )  of true path AB - start at E.
  const char* MethodName= "G4MultiLevelLocator::EstimateIntersectionPoint()";
  
  G4bool found_approximate_intersection = false;
  G4bool there_is_no_intersection       = false;

  G4FieldTrack  CurrentA_PointVelocity = CurveStartPointVelocity; 
  G4FieldTrack  CurrentB_PointVelocity = CurveEndPointVelocity;
  G4ThreeVector CurrentE_Point = TrialPoint;
  G4bool        validNormalAtE = false;
  G4ThreeVector NormalAtEntry;

  G4FieldTrack  ApproxIntersecPointV(CurveEndPointVelocity); // FT-Def-Construct
  G4bool        validIntersectP= true;   // Is it current ?
  G4double      NewSafety = 0.0;
  G4bool        last_AF_intersection = false;   

  auto    integrDriver = GetChordFinderFor()->GetIntegrationDriver();
  G4bool  driverReIntegrates = integrDriver->DoesReIntegrate();
  
  G4bool first_section = true;
  recalculatedEndPoint = false; 
  G4bool restoredFullEndpoint = false;

  unsigned int substep_no = 0;

  // Statistics for substeps
  static G4ThreadLocal unsigned int max_no_seen= 0;

#ifdef G4DEBUG_FIELD
  unsigned int trigger_substepno_print = 0;
  const G4double tolerance = 1.0e-8 * CLHEP::mm;
  unsigned int biggest_depth = 0;
    // using kInitialisingCL = G4LocatorChangeRecord::kInitialisingCL;
#endif

  // Log the location, iteration of changes in A,B 
  //----------------------------------------------
  static G4ThreadLocal G4LocatorChangeLogger endChangeA("StartPointA"),
    endChangeB("EndPointB"), recApproxPoint("ApproxPoint"),
    pointH_logger("Trial points 'E': position, normal");
  unsigned int eventCount = 0;

  if (fCheckMode)
  {
    // Clear previous call's data
    endChangeA.clear();
    endChangeB.clear();
    recApproxPoint.clear();
    pointH_logger.clear();

    // Record the initialisation
    ++eventCount;
    endChangeA.AddRecord( G4LocatorChangeRecord::kInitialisingCL, substep_no,
                          eventCount, CurrentA_PointVelocity );
    endChangeB.AddRecord( G4LocatorChangeRecord::kInitialisingCL, substep_no,
                          eventCount, CurrentB_PointVelocity );
  }

  //--------------------------------------------------------------------------  
  //  Algorithm for the case if progress in founding intersection is too slow.
  //  Process is defined too slow if after N=param_substeps advances on the
  //  path, it will be only 'fraction_done' of the total length.
  //  In this case the remaining length is divided in two half and 
  //  the loop is restarted for each half.  
  //  If progress is still too slow, the division in two halfs continue
  //  until 'max_depth'.
  //--------------------------------------------------------------------------

  const G4int param_substeps = 5;  // Test value for the maximum number
                                   // of substeps
  const G4double fraction_done = 0.3;

  G4bool Second_half = false;    // First half or second half of divided step

  // We need to know this for the 'final_section':
  // real 'final_section' or first half 'final_section'
  // In algorithm it is considered that the 'Second_half' is true
  // and it becomes false only if we are in the first-half of level
  // depthness or if we are in the first section

  unsigned int depth = 0; // Depth counts subdivisions of initial step made
  ++fNumCalls;
  
  NormalAtEntry = GetSurfaceNormal(CurrentE_Point, validNormalAtE);

  if (fCheckMode)
  {
    pointH_logger.AddRecord( G4LocatorChangeRecord::kInitialisingCL,
                             substep_no, eventCount,
                             G4FieldTrack( CurrentE_Point,0.,NormalAtEntry,0.,
                                           0., 1.,G4ThreeVector(0.),0.,0.) );
  #if (G4DEBUG_FIELD>1) 
    G4ThreeVector  StartPosition = CurveStartPointVelocity.GetPosition(); 
    if( (TrialPoint - StartPosition).mag2() < sqr(tolerance)) 
    {
       ReportImmediateHit( MethodName, StartPosition, TrialPoint,
                           tolerance, fNumCalls); 
    }
  #endif
  }

  // Intermediates Points on the Track = Subdivided Points must be stored.
  // Give the initial values to 'InterMedFt'
  // Important is 'ptrInterMedFT[0]', it saves the 'EndCurvePoint'
  //
  *ptrInterMedFT[0] = CurveEndPointVelocity;
  for ( auto idepth=1; idepth<max_depth+1; ++idepth )
  {
    *ptrInterMedFT[idepth] = CurveStartPointVelocity;
  }

  // Final_section boolean store
  //
  G4bool fin_section_depth[max_depth];
  for ( auto idepth=0; idepth<max_depth; ++idepth )
  {
    fin_section_depth[idepth] = true;
  }
  // 'SubStartPoint' is needed to calculate the length of the divided step
  //
  G4FieldTrack SubStart_PointVelocity = CurveStartPointVelocity;
   
  do  // Loop checking, 07.10.2016, J.Apostolakis
  {
    unsigned int substep_no_p = 0;
    G4bool sub_final_section = false; // the same as final_section,
                                      // but for 'sub_section'
    SubStart_PointVelocity = CurrentA_PointVelocity;
 
    do // Loop checking, 07.10.2016, J.Apostolakis
    { // REPEAT param
      G4ThreeVector Point_A = CurrentA_PointVelocity.GetPosition();  
      G4ThreeVector Point_B = CurrentB_PointVelocity.GetPosition();
       
#ifdef G4DEBUG_FIELD
      const G4double lenA = CurrentA_PointVelocity.GetCurveLength() ;
      const G4double lenB = CurrentB_PointVelocity.GetCurveLength() ;            
      G4double curv_lenAB = lenB - lenA;
      G4double     distAB = (Point_B - Point_A).mag();
      if( curv_lenAB < distAB * ( 1. - 10.*fiEpsilonStep ) )
      {
        G4cerr << "ERROR> (Start) Point A coincides with or has gone past (end) point B"           
               << "MLL: iters = " << substep_no << G4endl;
        G4long op=G4cerr.precision(6);
        G4cerr << "       Difference = " << distAB - curv_lenAB
               << " exceeds limit of relative dist (10*epsilon)= " << 10*fiEpsilonStep
               << "  i.e. limit = " << 10 * fiEpsilonStep * distAB << G4endl;
        G4cerr.precision(9);        
        G4cerr << "        Len A, B = " << lenA << " " << lenB << G4endl
               << "        Position A: " << Point_A << G4endl
               << "        Position B: " << Point_B << G4endl;
        G4cerr.precision(op);
        // G4LocatorChangeRecord::ReportVector(G4cerr, "endPointB", endChangeB );
        // G4cerr<<"EndPoints A(start) and B(end): combined changes " << G4endl;
        if (fCheckMode) {
           G4LocatorChangeLogger::ReportEndChanges(G4cerr, endChangeA, endChangeB);
        }
      }
#endif    
      if( !validIntersectP ){
        G4ExceptionDescription errmsg;
        errmsg << "Assertion FAILURE - invalid (stale) Interection point. Substep: "
               << substep_no << " call: " << fNumCalls << G4endl;
        if (fCheckMode)
           G4LocatorChangeRecord::ReportEndChanges(errmsg, endChangeA, endChangeB );
        G4Exception("G4MultiLevelLocator::EstimateIntersectionPoint", "GeomNav0004",
                    JustWarning, errmsg);
      }
      
      // F = a point on true AB path close to point E 
      // (the closest if possible)
      //
      ApproxIntersecPointV = GetChordFinderFor()
                             ->ApproxCurvePointV( CurrentA_PointVelocity, 
                                                  CurrentB_PointVelocity, 
                                                  CurrentE_Point,
                                                  GetEpsilonStepFor() );
        // The above method is the key & most intuitive part ...

#ifdef G4DEBUG_FIELD      
      recApproxPoint.push_back(G4LocatorChangeRecord(G4LocatorChangeRecord::kInvalidCL,
                                                  substep_no, eventCount, ApproxIntersecPointV ) );
      G4double lenIntsc= ApproxIntersecPointV.GetCurveLength();
      G4double checkVsEnd= lenB - lenIntsc;

      if( lenIntsc > lenB ) 
      {
        std::ostringstream errmsg;
        errmsg.precision(17);
        G4double ratio    = checkVsEnd / lenB;
        G4double ratioTol = std::fabs(ratio) / tolerance;
        errmsg << "Intermediate F point is past end B point" << G4endl
        << "   l( intersection ) = " << lenIntsc << G4endl
        << "   l( endpoint     ) = " << lenB     << G4endl;
        errmsg.precision(8);
        errmsg << "   l_end - l_inters  = " << checkVsEnd << G4endl
        << "          / l_end      = " << ratio << G4endl
        << "   ratio  / tolerance  = " << ratioTol  << G4endl;
        if( ratioTol < 1.0 )
          G4Exception(MethodName, "GeomNav0003", JustWarning, errmsg );
        else
          G4Exception(MethodName, "GeomNav0003", FatalException, errmsg );
      }
#endif

      G4ThreeVector CurrentF_Point= ApproxIntersecPointV.GetPosition();

      // First check whether EF is small - then F is a good approx. point 
      // Calculate the length and direction of the chord AF
      //
      G4ThreeVector  ChordEF_Vector = CurrentF_Point - CurrentE_Point;

      G4ThreeVector  NewMomentumDir = ApproxIntersecPointV.GetMomentumDir(); 
      G4double       MomDir_dot_Norm = NewMomentumDir.dot( NormalAtEntry );

#ifdef G4DEBUG_FIELD
      if( fVerboseLevel > 3 )
      { 
         G4ThreeVector  ChordAB           = Point_B - Point_A;
         G4double       ABchord_length    = ChordAB.mag(); 
         G4double       MomDir_dot_ABchord;
         MomDir_dot_ABchord = (1.0 / ABchord_length)
                            * NewMomentumDir.dot( ChordAB );
         G4VIntersectionLocator::ReportTrialStep( substep_no, ChordAB,
              ChordEF_Vector, NewMomentumDir, NormalAtEntry, validNormalAtE ); 
         G4cout << " dot( MomentumDir, ABchord_unit ) = "
                << MomDir_dot_ABchord << G4endl;
      }
#endif
      G4bool adequate_angle =
             ( MomDir_dot_Norm >= 0.0 ) // Can use ( > -epsilon) instead
          || (! validNormalAtE );       // Invalid, cannot use this criterion
      G4double EF_dist2 = ChordEF_Vector.mag2();
      if ( ( EF_dist2 <= sqr(fiDeltaIntersection) && ( adequate_angle ) )
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
          G4ThreeVector MomentumDir=ApproxIntersecPointV.GetMomentumDirection();
          G4bool goodCorrection = AdjustmentOfFoundIntersection(Point_A,
                                    CurrentE_Point, CurrentF_Point, MomentumDir,
                                    last_AF_intersection, IP, NewSafety,
                                    previousSafety, previousSftOrigin );
          if ( goodCorrection )
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
        G4bool Intersects_FB = false;
        G4bool Intersects_AF = IntersectChord( Point_A,   CurrentF_Point,
                                               NewSafety, previousSafety,
                                               previousSftOrigin,
                                               stepLengthAF,
                                               PointG );
        last_AF_intersection = Intersects_AF;
        if( Intersects_AF )
        {
          // G is our new Candidate for the intersection point.
          // It replaces  "E" and we will see if it's good enough.
          CurrentB_PointVelocity = ApproxIntersecPointV;  // B  <- F
          CurrentE_Point = PointG;                        // E  <- G

          validIntersectP = true;    // 'E' has been updated.

          G4bool validNormalLast; 
          NormalAtEntry  = GetSurfaceNormal( PointG, validNormalLast ); 
          validNormalAtE = validNormalLast; 

          // As we move point B, must take care in case the current
          // AF has no intersection to try current FB!!
          fin_section_depth[depth] = false;

          if (fCheckMode)
          {
            ++eventCount;
            endChangeB.push_back(
              G4LocatorChangeRecord(G4LocatorChangeRecord::kIntersectsAF,
                             substep_no, eventCount, CurrentB_PointVelocity) );
          }
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

          // Check whether any volumes are encountered by the chord FB
          // ---------------------------------------------------------

          Intersects_FB = IntersectChord( CurrentF_Point, Point_B, 
                                          NewSafety, previousSafety,
                                          previousSftOrigin,
                                          stepLengthFB,
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
            //
            CurrentA_PointVelocity = ApproxIntersecPointV;
            CurrentE_Point = PointH;

            validIntersectP = true;    // 'E' has been updated.

            G4bool validNormalH;
            NormalAtEntry  = GetSurfaceNormal( PointH, validNormalH ); 
            validNormalAtE = validNormalH;

            if (fCheckMode)
            {
              ++eventCount;
              endChangeA.push_back( 
                 G4LocatorChangeRecord(G4LocatorChangeRecord::kIntersectsFB,
                             substep_no, eventCount, CurrentA_PointVelocity) );
              G4FieldTrack intersectH_pn('0');  // Point and normal
                                                // nothing else will be valid
              intersectH_pn.SetPosition( PointH );
              intersectH_pn.SetMomentum( NormalAtEntry ); 
              pointH_logger.AddRecord(G4LocatorChangeRecord::kIntersectsFB,
                                    substep_no, eventCount, intersectH_pn );
            }
          }
          else  // not Intersects_FB
          {
            validIntersectP = false;    // Intersections are now stale
            if( fin_section_depth[depth] )
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
                Second_half = true; 
                sub_final_section = true;
              }
            }
            else
            {
              CurrentA_PointVelocity = CurrentB_PointVelocity;  // Got to B
              CurrentB_PointVelocity = (depth==0) ? CurveEndPointVelocity 
                                                  : *ptrInterMedFT[depth] ;
              SubStart_PointVelocity = CurrentA_PointVelocity;
              restoredFullEndpoint = true;

              validIntersectP = false;  // 'E' has NOT been updated.

              if (fCheckMode)
              {
                ++eventCount;
                endChangeA.push_back(
                  G4LocatorChangeRecord(
                    G4LocatorChangeRecord::kNoIntersectAForFB,
                             substep_no, eventCount, CurrentA_PointVelocity) ); 
                endChangeB.push_back(
                  G4LocatorChangeRecord (
                    G4LocatorChangeRecord::kNoIntersectAForFB,
                             substep_no, eventCount, CurrentB_PointVelocity) );
              }
            }
          } // Endif (Intersects_FB)
        } // Endif (Intersects_AF)

        G4int errorEndPt = 0; // Default: no error (if not calling CheckAnd...

        G4bool recalculatedB= false;
        if( driverReIntegrates )
        {
          G4FieldTrack RevisedB_FT = CurrentB_PointVelocity;
          recalculatedB= CheckAndReEstimateEndpoint(CurrentA_PointVelocity,
                                                    CurrentB_PointVelocity,
                                                    RevisedB_FT,
                                                    errorEndPt );
          if( recalculatedB )
          {
            CurrentB_PointVelocity = RevisedB_FT;  // Use it !
            // Do not invalidate intersection F -- it is still roughly OK.
            //
            // The best course would be to invalidate (reset validIntersectP)
            // BUT if we invalidate it, we must re-estimate it somewhere! E.g.
            //   validIntersectP = false;  // 'E' has NOT been updated.
 
            if ( (fin_section_depth[depth])           // real final section
               &&( first_section || ((Second_half)&&(depth==0)) ) )
            {
              recalculatedEndPoint = true;
              IntersectedOrRecalculatedFT = RevisedB_FT;
              // So that we can return it, if it is the endpoint!
            }
            // else
            //  Move forward the other points 
            //   - or better flag it, so that they are re-computed when next used
            //     [ Implementation: a counter for # of recomputations
            //       => avoids extra work]
          }
          if (fCheckMode)
          {
            ++eventCount;
            endChangeB.push_back(
              G4LocatorChangeRecord( G4LocatorChangeRecord::kRecalculatedB,
                                     substep_no, eventCount, RevisedB_FT ) );
          }
        }
        else
        {
          if( CurrentB_PointVelocity.GetCurveLength() < CurrentA_PointVelocity.GetCurveLength() )
            errorEndPt = 2;
        }

        if( errorEndPt > 1 )  // errorEndPt = 1 is milder, just: len(B)=len(A)
        {
          std::ostringstream errmsg;

          ReportReversedPoints(errmsg,
                               CurveStartPointVelocity, CurveEndPointVelocity,
                               NewSafety, fiEpsilonStep,
                               CurrentA_PointVelocity, CurrentB_PointVelocity,
                               SubStart_PointVelocity, CurrentE_Point,
                               ApproxIntersecPointV, substep_no, substep_no_p, depth);
          
          if (fCheckMode) {
            G4LocatorChangeRecord::ReportEndChanges(errmsg, endChangeA, endChangeB );
          }

          errmsg << G4endl << " * Location: " << MethodName
                 << "- After EndIf(Intersects_AF)" << G4endl;
          errmsg << " * Bool flags:  Recalculated = " << recalculatedB
                 << "   Intersects_AF = " << Intersects_AF
                 << "   Intersects_FB = " << Intersects_FB << G4endl;
          errmsg << " * Number of calls to MLL:EIP= " << fNumCalls << G4endl;
          G4Exception(MethodName, "GeomNav0003", FatalException, errmsg);
        }
        if( restoredFullEndpoint )
        {
          fin_section_depth[depth] = restoredFullEndpoint;
          restoredFullEndpoint = false;
        }
      } // EndIf ( E is close enough to the curve, ie point F. )
        // tests ChordAF_Vector.mag() <= maximum_lateral_displacement 

#ifdef G4DEBUG_FIELD
      if( trigger_substepno_print == 0)
      {
        trigger_substepno_print= fWarnSteps - 20;
      }

      if( substep_no >= trigger_substepno_print )
      {
        G4cout << "Difficulty in converging in " << MethodName
               << "    Substep no = " << substep_no << G4endl;
        if( substep_no == trigger_substepno_print )
        {
          G4cout << " Start:   ";
          printStatus( CurveStartPointVelocity, CurveEndPointVelocity,
                       -1.0, NewSafety, 0 );
          if( fCheckMode ) {
            G4LocatorChangeRecord::ReportEndChanges(G4cout, endChangeA, endChangeB );
          } else {
            G4cout << " ** For more information enable 'check mode' in G4MultiLevelLocator "
                   << "-- (it saves and can output change records) " << G4endl;
          }
        }
        G4cout << " Point A: "; 
        printStatus( CurrentA_PointVelocity, CurrentA_PointVelocity,
                     -1.0, NewSafety, substep_no-1 );
        G4cout << " Point B: ";
        printStatus( CurrentA_PointVelocity, CurrentB_PointVelocity,
                     -1.0, NewSafety, substep_no );
      }
#endif
      ++substep_no; 
      ++substep_no_p;

    } while (  ( ! found_approximate_intersection )
            && ( ! there_is_no_intersection )     
            && validIntersectP        // New condition:  must refresh intersection !!
            && ( substep_no_p <= param_substeps) );  // UNTIL found or
                                                     // failed param substep

    if( (!found_approximate_intersection) && (!there_is_no_intersection) )
    {
      G4double did_len = std::abs( CurrentA_PointVelocity.GetCurveLength()
                       - SubStart_PointVelocity.GetCurveLength()); 
      G4double all_len = std::abs( CurrentB_PointVelocity.GetCurveLength()
                       - SubStart_PointVelocity.GetCurveLength());
   
      G4double distAB = -1;
      //
      // Is progress is too slow, and is it possible to go deeper?
      // If so, then *halve the step*
      //              ==============
      if(   (did_len < fraction_done*all_len)
         && (depth<max_depth) && (!sub_final_section) )
      {
#ifdef G4DEBUG_FIELD         
        static G4ThreadLocal unsigned int numSplits=0;   // For debugging only 
        biggest_depth = std::max(depth, biggest_depth);
        ++numSplits;
#endif
        Second_half = false;
        ++depth;
        first_section = false;

        G4double Sub_len = (all_len-did_len)/(2.);
        G4FieldTrack midPoint = CurrentA_PointVelocity;
        G4bool fullAdvance=              
           integrDriver->AccurateAdvance(midPoint, Sub_len, fiEpsilonStep);
                         
        ++fNumAdvanceTrials;
        if( fullAdvance )  { ++fNumAdvanceFull; }

        G4double lenAchieved=
           midPoint.GetCurveLength()-CurrentA_PointVelocity.GetCurveLength();
        
        const G4double adequateFraction = (1.0-CLHEP::perThousand);
        G4bool goodAdvance = (lenAchieved >= adequateFraction * Sub_len);
        if ( goodAdvance )  { ++fNumAdvanceGood; }

#ifdef G4DEBUG_FIELD
        else  //  !goodAdvance
        {
           G4cout << "MLL> AdvanceChordLimited not full at depth=" << depth 
                  << "  total length achieved = " << lenAchieved << " of "
                  << Sub_len << "  fraction= ";
           if (Sub_len != 0.0 ) { G4cout << lenAchieved / Sub_len; }
           else                 { G4cout << "DivByZero"; }
           G4cout << "  Good-enough fraction = " << adequateFraction;
           G4cout << "  Number of call to mll = " << fNumCalls
                  << "   iteration " << substep_no
                  << "  inner =  " << substep_no_p << G4endl;
           G4cout << "  Epsilon of integration = " << fiEpsilonStep << G4endl;
           G4cout << "  State at start is      = " << CurrentA_PointVelocity
                  << G4endl
                  << "        at end (midpoint)= " << midPoint << G4endl;
           G4cout << "  Particle mass = " << midPoint.GetRestMass() << G4endl;

           G4EquationOfMotion *equation = integrDriver->GetEquationOfMotion();
           ReportFieldValue( CurrentA_PointVelocity, "start", equation );
           ReportFieldValue( midPoint, "midPoint", equation );            
           G4cout << "  Original Start = "
                  << CurveStartPointVelocity << G4endl;
           G4cout << "  Original   End = "
                  << CurveEndPointVelocity << G4endl;
           G4cout << "  Original TrialPoint = "
                  << TrialPoint << G4endl;
           G4cout << "  (this is STRICT mode) "
                  << "  num Splits= " << numSplits;           
           G4cout << G4endl;
        }
#endif

        *ptrInterMedFT[depth] =  midPoint;
        CurrentB_PointVelocity = midPoint; 

        if (fCheckMode)
        {
          ++eventCount;
          endChangeB.push_back(
            G4LocatorChangeRecord( G4LocatorChangeRecord::kInsertingMidPoint,
                                   substep_no, eventCount, midPoint) );
        }        
        
        // Adjust 'SubStartPoint' to calculate the 'did_length' in next loop
        //
        SubStart_PointVelocity = CurrentA_PointVelocity;

        // Find new trial intersection point needed at start of the loop
        //
        G4ThreeVector Point_A = CurrentA_PointVelocity.GetPosition();
        G4ThreeVector Point_B = CurrentB_PointVelocity.GetPosition();   

        G4ThreeVector PointGe;
        GetNavigatorFor()->LocateGlobalPointWithinVolume(Point_A);
        G4bool Intersects_AB = IntersectChord(Point_A, Point_B,
                                              NewSafety, previousSafety,
                                              previousSftOrigin, distAB,
                                              PointGe);
        if( Intersects_AB )
        {
          last_AF_intersection = Intersects_AB;
          CurrentE_Point = PointGe;
          fin_section_depth[depth] = true;

          validIntersectP = true;  // 'E' has been updated.

          G4bool validNormalAB; 
          NormalAtEntry  = GetSurfaceNormal( PointGe, validNormalAB ); 
          validNormalAtE = validNormalAB;
        }
        else
        {
          // No intersection found for first part of curve
          // (CurrentA,InterMedPoint[depth]). Go to the second part
          //
          Second_half = true;

          validIntersectP=  false;  // No new 'E' chord intersection found
        }
      } // if did_len

      G4bool unfinished = Second_half;
      while ( unfinished && (depth>0) )  // Loop checking, 07.10.2016, JA
      {
        // Second part of curve (InterMed[depth],Intermed[depth-1])) 
        // On the depth-1 level normally we are on the 'second_half'

        //  Find new trial intersection point needed at start of the loop
        //
        SubStart_PointVelocity = *ptrInterMedFT[depth];
        CurrentA_PointVelocity = *ptrInterMedFT[depth];
        CurrentB_PointVelocity = *ptrInterMedFT[depth-1];

        if (fCheckMode)
        {
          ++eventCount;
          G4LocatorChangeRecord chngPop_a( G4LocatorChangeRecord::kLevelPop,
                              substep_no, eventCount, CurrentA_PointVelocity);
          endChangeA.push_back( chngPop_a );
          G4LocatorChangeRecord chngPop_b( G4LocatorChangeRecord::kLevelPop,
                              substep_no, eventCount, CurrentB_PointVelocity);
          endChangeB.push_back( chngPop_b );
        }
        
        // Ensure that the new endpoints are not further apart in space
        // than on the curve due to different errors in the integration
        //
        G4int  errorEndPt = -1;
        G4bool recalculatedB= false;
        if( driverReIntegrates )
        {
          // Ensure that the new endpoints are not further apart in space
          // than on the curve due to different errors in the integration
          //
          G4FieldTrack RevisedEndPointFT = CurrentB_PointVelocity;
          recalculatedB = 
             CheckAndReEstimateEndpoint( CurrentA_PointVelocity,  
                                         CurrentB_PointVelocity,
                                         RevisedEndPointFT,
                                         errorEndPt );
          if( recalculatedB )
          {
            CurrentB_PointVelocity = RevisedEndPointFT;  // Use it !
            
            if ( depth == 1 )
            {
               recalculatedEndPoint = true;
               IntersectedOrRecalculatedFT = RevisedEndPointFT;
               // So that we can return it, if it is the endpoint!
            }
          }
          else
          {          
            if( CurrentB_PointVelocity.GetCurveLength() < CurrentA_PointVelocity.GetCurveLength() )
              errorEndPt = 2;
          }
          
          if (fCheckMode)
          {
            ++eventCount;
            endChangeB.push_back(
              G4LocatorChangeRecord(G4LocatorChangeRecord::kRecalculatedBagn,
                                    substep_no, eventCount, RevisedEndPointFT));
          }
        }
        if( errorEndPt > 1 )  // errorEndPt = 1 is milder, just: len(B)=len(A)
        {
          std::ostringstream errmsg; 
          ReportReversedPoints(errmsg, 
                    CurveStartPointVelocity, CurveEndPointVelocity,
                    NewSafety, fiEpsilonStep, 
                    CurrentA_PointVelocity, CurrentB_PointVelocity,
                    SubStart_PointVelocity, CurrentE_Point,
                    ApproxIntersecPointV, substep_no, substep_no_p, depth);
          errmsg << " * Location:  " << MethodName << "- Second-Half" << G4endl;
          errmsg << " * Recalculated = " << recalculatedB << G4endl; // false
          G4Exception(MethodName, "GeomNav0003", FatalException, errmsg);
        }
        
        G4ThreeVector Point_A = CurrentA_PointVelocity.GetPosition();
        G4ThreeVector Point_B = CurrentB_PointVelocity.GetPosition();
        G4ThreeVector PointGi;     
        GetNavigatorFor()->LocateGlobalPointWithinVolume(Point_A);
        G4bool Intersects_AB = IntersectChord(Point_A, Point_B, NewSafety,
                                              previousSafety,
                                              previousSftOrigin, distAB,
                                              PointGi);
        if( Intersects_AB )
        {
          last_AF_intersection = Intersects_AB;
          CurrentE_Point = PointGi;

          validIntersectP = true;  // 'E' has been updated.
          NormalAtEntry  = GetSurfaceNormal( PointGi, validNormalAtE ); 
        }
        else
        {
          validIntersectP = false;  // No new 'E' chord intersection found
          if( depth == 1)
          {
            there_is_no_intersection = true;
          }
        }
        depth--;
        fin_section_depth[depth] = true;
        unfinished = !validIntersectP;
      }
#ifdef G4DEBUG_FIELD
      if( ! ( validIntersectP || there_is_no_intersection ) )
      {
         // What happened ??
         G4cout << "MLL - WARNING Potential FAILURE: Conditions not met!"
                << G4endl
                << " Depth = " << depth << G4endl
                << " Num Substeps= " << substep_no << G4endl;
         G4cout << " Found intersection= " << found_approximate_intersection
                << G4endl;
         G4cout << " Progress report: -- " << G4endl;
         ReportProgress(G4cout, 
                        CurveStartPointVelocity, CurveEndPointVelocity,
                        substep_no, CurrentA_PointVelocity,
                        CurrentB_PointVelocity,
                        NewSafety, depth ); 
         G4cout << G4endl; 
      }
#endif      
    }  // if(!found_aproximate_intersection)

    assert( validIntersectP || there_is_no_intersection
                            || found_approximate_intersection);

  } while ( ( ! found_approximate_intersection )
            && ( ! there_is_no_intersection )     
            && ( substep_no <= fMaxSteps) ); // UNTIL found or failed

  if( substep_no > max_no_seen )
  {
    max_no_seen = substep_no; 
#ifdef G4DEBUG_FIELD
    if( max_no_seen > fWarnSteps )
    {
      trigger_substepno_print = max_no_seen-20; // Want to see last 20 steps 
    } 
#endif
  }

  if( !there_is_no_intersection && !found_approximate_intersection )
  {
     if( substep_no >= fMaxSteps)
     {
        // Since we cannot go further (yet), we return as far as we have gone

        recalculatedEndPoint = true;
        IntersectedOrRecalculatedFT = CurrentA_PointVelocity;
        found_approximate_intersection = false; 
     
        std::ostringstream message;
        message << G4endl;
        message << "Convergence is requiring too many substeps: "
                << substep_no << "  ( limit = "<<  fMaxSteps << ")"
                << G4endl 
                << "  Abandoning effort to intersect. " << G4endl << G4endl;
        message << "    Number of calls to MLL: " <<   fNumCalls;
        message << "  iteration = " << substep_no <<G4endl << G4endl;
        
        message.precision( 12 ); 
        G4double done_len = CurrentA_PointVelocity.GetCurveLength(); 
        G4double full_len = CurveEndPointVelocity.GetCurveLength();
        message << "        Undertaken only length: " << done_len
                << " out of " << full_len << " required." << G4endl
                << "        Remaining length = " << full_len - done_len; 
        
        message << "     Start and end-point of requested Step:" << G4endl;
        printStatus( CurveStartPointVelocity, CurveEndPointVelocity,
                     -1.0, NewSafety, 0,     message, -1 );
        message << "        Start and end-point of current Sub-Step:" << G4endl;
        printStatus( CurrentA_PointVelocity, CurrentA_PointVelocity,
                     -1.0, NewSafety, substep_no-1, message, -1 );
        printStatus( CurrentA_PointVelocity, CurrentB_PointVelocity,
                     -1.0, NewSafety, substep_no, message, -1 );
        
        G4Exception(MethodName, "GeomNav0003", JustWarning, message);
     }
     else if( substep_no >= fWarnSteps)
     {  
        std::ostringstream message;
        message << "Many substeps while trying to locate intersection."
                << G4endl
                << "          Undertaken length: "  
                << CurrentB_PointVelocity.GetCurveLength()
                << " - Needed: "  << substep_no << " substeps." << G4endl
                << "          Warning number = " << fWarnSteps
                << " and maximum substeps = " << fMaxSteps;
        G4Exception(MethodName, "GeomNav1002", JustWarning, message);
     }
  }
  
  return  (!there_is_no_intersection) && found_approximate_intersection;
    //  Success or failure
}

void G4MultiLevelLocator::ReportStatistics()
{
   G4cout << " Number of calls = " << fNumCalls << G4endl;
   G4cout << " Number of split level ('advances'):  "
          << fNumAdvanceTrials << G4endl;
   G4cout << " Number of full advances:             "
          << fNumAdvanceGood << G4endl;
   G4cout << " Number of good advances:             "
          << fNumAdvanceFull << G4endl;
}

void G4MultiLevelLocator::ReportFieldValue( const G4FieldTrack& locationPV,
                                            const char* nameLoc,
                                            const G4EquationOfMotion* equation )
{
   enum { maxNumFieldComp = 24 };

   G4ThreeVector position = locationPV.GetPosition();   
   G4double startPoint[4] = { position.x(), position.y(), position.z(),
                              locationPV.GetLabTimeOfFlight() };
   G4double FieldVec[maxNumFieldComp]; // 24 ;
   for (auto i=0; i<maxNumFieldComp; ++i )
   {
     FieldVec[i] = 0.0;
   }
   equation->GetFieldValue( startPoint, FieldVec);
   G4cout << "  B-field value (" << nameLoc << ")=   "
          << FieldVec[0] << " " << FieldVec[1] << " " << FieldVec[2];
   G4double Emag2= G4ThreeVector( FieldVec[3],
                                  FieldVec[4],
                                  FieldVec[5] ).mag2(); 
   if( Emag2 > 0.0 )
   {
      G4cout << " Electric = " << FieldVec[3] << " "
                               << FieldVec[4] << " "
                               << FieldVec[5]<< G4endl;
   }
   return;
}
