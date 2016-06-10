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
// $Id: G4ChordFinder.cc 95290 2016-02-04 08:14:10Z gcosmo $
//
//
// 25.02.97 - John Apostolakis - Design and implementation 
// -------------------------------------------------------------------

#include <iomanip>

#include "G4ChordFinder.hh"
#include "G4SystemOfUnits.hh"
#include "G4MagneticField.hh"
#include "G4Mag_UsualEqRhs.hh"
#include "G4ClassicalRK4.hh"


// ..........................................................................

G4ChordFinder::G4ChordFinder(G4MagInt_Driver* pIntegrationDriver)
  : fDefaultDeltaChord( 0.25 * mm ),      // Parameters
    fDeltaChord( fDefaultDeltaChord ),    //   Internal parameters
    fFirstFraction(0.999), fFractionLast(1.00),  fFractionNextEstimate(0.98), 
    fMultipleRadius(15.0), 
    fStatsVerbose(0),
    fDriversStepper(0),                    // Dependent objects 
    fAllocatedStepper(false),
    fEquation(0),      
    fTotalNoTrials_FNC(0), fNoCalls_FNC(0), fmaxTrials_FNC(0)
{
  // Simple constructor -- it does not create equation
  fIntgrDriver= pIntegrationDriver;
  fAllocatedStepper= false;

  fLastStepEstimate_Unconstrained = DBL_MAX;          // Should move q, p to

  SetFractions_Last_Next( fFractionLast, fFractionNextEstimate);  
    // check the values and set the other parameters
}


// ..........................................................................

G4ChordFinder::G4ChordFinder( G4MagneticField*        theMagField,
                              G4double                stepMinimum, 
                              G4MagIntegratorStepper* pItsStepper )
  : fDefaultDeltaChord( 0.25 * mm ),     // Constants 
    fDeltaChord( fDefaultDeltaChord ),   // Parameters
    fFirstFraction(0.999), fFractionLast(1.00),  fFractionNextEstimate(0.98), 
    fMultipleRadius(15.0), 
    fStatsVerbose(0),
    fDriversStepper(0),                  //  Dependent objects
    fAllocatedStepper(false),
    fEquation(0), 
    fTotalNoTrials_FNC(0), fNoCalls_FNC(0), fmaxTrials_FNC(0)  // State - stats
{
  //  Construct the Chord Finder
  //  by creating in inverse order the  Driver, the Stepper and EqRhs ...

  G4Mag_EqRhs *pEquation = new G4Mag_UsualEqRhs(theMagField);
  fEquation = pEquation;                            
  fLastStepEstimate_Unconstrained = DBL_MAX;          // Should move q, p to
                                                     //    G4FieldTrack ??

  SetFractions_Last_Next( fFractionLast, fFractionNextEstimate);  
    // check the values and set the other parameters

  // --->>  Charge    Q = 0 
  // --->>  Momentum  P = 1       NOMINAL VALUES !!!!!!!!!!!!!!!!!!

  if( pItsStepper == 0 )
  { 
     pItsStepper = fDriversStepper = new G4ClassicalRK4(pEquation);
     fAllocatedStepper= true;
  }
  else
  {
     fAllocatedStepper= false; 
  }
  fIntgrDriver = new G4MagInt_Driver(stepMinimum, pItsStepper, 
                                     pItsStepper->GetNumberOfVariables() );
}


// ......................................................................

G4ChordFinder::~G4ChordFinder()
{
  delete   fEquation; // fIntgrDriver->pIntStepper->theEquation_Rhs;
  if( fAllocatedStepper)
  { 
     delete fDriversStepper; 
  }
  delete   fIntgrDriver; 

  if( fStatsVerbose ) { PrintStatistics(); }
}


// ......................................................................

void   
G4ChordFinder::SetFractions_Last_Next( G4double fractLast, G4double fractNext )
{ 
  // Use -1.0 as request for Default.
  if( fractLast == -1.0 )   fractLast = 1.0;   // 0.9;
  if( fractNext == -1.0 )   fractNext = 0.98;  // 0.9; 

  // fFirstFraction  = 0.999; // Orig 0.999 A safe value, range: ~ 0.95 - 0.999
  // fMultipleRadius = 15.0;  // For later use, range: ~  2 - 20 

  if( fStatsVerbose )
  { 
    G4cout << " ChordFnd> Trying to set fractions: "
           << " first " << fFirstFraction
           << " last " <<  fractLast
           << " next " <<  fractNext
           << " and multiple " << fMultipleRadius
           << G4endl;
  } 

  if( (fractLast > 0.0) && (fractLast <=1.0) ) 
  {
    fFractionLast= fractLast;
  }
  else
  {
    G4cerr << "G4ChordFinder::SetFractions_Last_Next: Invalid "
           << " fraction Last = " << fractLast
           << " must be  0 <  fractionLast <= 1 " << G4endl;
  }
  if( (fractNext > 0.0) && (fractNext <1.0) )
  {
    fFractionNextEstimate = fractNext;
  }
  else
  {
    G4cerr << "G4ChordFinder:: SetFractions_Last_Next: Invalid "
           << " fraction Next = " << fractNext
           << " must be  0 <  fractionNext < 1 " << G4endl;
  }
}


// ......................................................................

G4double 
G4ChordFinder::AdvanceChordLimited( G4FieldTrack& yCurrent,
                                    G4double      stepMax,
                                    G4double      epsStep,
                                    const G4ThreeVector latestSafetyOrigin,
                                    G4double       latestSafetyRadius )
{
  G4double stepPossible;
  G4double dyErr;
  G4FieldTrack yEnd( yCurrent);
  G4double  startCurveLen= yCurrent.GetCurveLength();
  G4double nextStep;
  //            *************
  stepPossible= FindNextChord(yCurrent, stepMax, yEnd, dyErr, epsStep,
                              &nextStep, latestSafetyOrigin, latestSafetyRadius
                             );
  //            *************

  G4bool good_advance;

  if ( dyErr < epsStep * stepPossible )
  {
     // Accept this accuracy.

     yCurrent = yEnd;
     good_advance = true; 
  }
  else
  {  
     // Advance more accurately to "end of chord"
     //                           ***************
     good_advance = fIntgrDriver->AccurateAdvance(yCurrent, stepPossible,
                                                  epsStep, nextStep);
     if ( ! good_advance )
     { 
       // In this case the driver could not do the full distance
       stepPossible= yCurrent.GetCurveLength()-startCurveLen;
     }
  }
  return stepPossible;
}


// ............................................................................

G4double
G4ChordFinder::FindNextChord( const  G4FieldTrack& yStart,
                                     G4double     stepMax,
                                     G4FieldTrack&   yEnd, // Endpoint
                                     G4double&   dyErrPos, // Error of endpoint
                                     G4double    epsStep,
                                     G4double*  pStepForAccuracy, 
                              const  G4ThreeVector, //  latestSafetyOrigin,
                                     G4double       //  latestSafetyRadius 
                                        )
{
  // Returns Length of Step taken

  G4FieldTrack yCurrent=  yStart;  
  G4double    stepTrial, stepForAccuracy;
  G4double    dydx[G4FieldTrack::ncompSVEC]; 

  //  1.)  Try to "leap" to end of interval
  //  2.)  Evaluate if resulting chord gives d_chord that is good enough.
  // 2a.)  If d_chord is not good enough, find one that is.
  
  G4bool    validEndPoint= false;
  G4double  dChordStep, lastStepLength; //  stepOfLastGoodChord;

  fIntgrDriver-> GetDerivatives( yCurrent, dydx );

  unsigned int        noTrials=0;
  const unsigned int  maxTrials= 300; // Avoid endless loop for bad convergence 

  const G4double safetyFactor= fFirstFraction; //  0.975 or 0.99 ? was 0.999

  stepTrial = std::min( stepMax, safetyFactor*fLastStepEstimate_Unconstrained );

  G4double newStepEst_Uncons= 0.0; 
  G4double stepForChord;
  do
  { 
     yCurrent = yStart;    // Always start from initial point
    
     //            ************
     fIntgrDriver->QuickAdvance( yCurrent, dydx, stepTrial, 
                                 dChordStep, dyErrPos);
     //            ************
     
     //  We check whether the criterion is met here.
     validEndPoint = AcceptableMissDist(dChordStep);

     lastStepLength = stepTrial; 

     // This method estimates to step size for a good chord.
     stepForChord = NewStep(stepTrial, dChordStep, newStepEst_Uncons );

     if( ! validEndPoint )
     {
        if( stepTrial<=0.0 )
        {
          stepTrial = stepForChord;
        }
        else if (stepForChord <= stepTrial)
        {
          // Reduce by a fraction, possibly up to 20% 
          stepTrial = std::min( stepForChord, fFractionLast * stepTrial);
        }
        else
        {
          stepTrial *= 0.1;
        }
     }
     noTrials++; 
  }
  while( (! validEndPoint) && (noTrials < maxTrials) );   // End of do-while  RKD 

  if( noTrials >= maxTrials )
  {
      std::ostringstream message;
      message << "Exceeded maximum number of trials= " << maxTrials << G4endl
              << "Current sagita dist= " << dChordStep << G4endl
              << "Step sizes (actual and proposed): " << G4endl
              << "Last trial =         " << lastStepLength  << G4endl
              << "Next trial =         " << stepTrial  << G4endl
              << "Proposed for chord = " << stepForChord  << G4endl              
              ;
      G4Exception("G4ChordFinder::FindNextChord()", "GeomField0003",
                  JustWarning, message);
  }

  if( newStepEst_Uncons > 0.0  )
  {
     fLastStepEstimate_Unconstrained= newStepEst_Uncons;
  }

  AccumulateStatistics( noTrials );

  if( pStepForAccuracy )
  { 
     // Calculate the step size required for accuracy, if it is needed
     //
     G4double dyErr_relative = dyErrPos/(epsStep*lastStepLength);
     if( dyErr_relative > 1.0 )
     {
        stepForAccuracy = fIntgrDriver->ComputeNewStepSize( dyErr_relative,
                                                            lastStepLength );
     }
     else
     {
        stepForAccuracy = 0.0;   // Convention to show step was ok 
     }
     *pStepForAccuracy = stepForAccuracy;
  }

#ifdef  TEST_CHORD_PRINT
  static int dbg=0;
  if( dbg )
  {
    G4cout << "ChordF/FindNextChord:  NoTrials= " << noTrials 
           << " StepForGoodChord=" << std::setw(10) << stepTrial << G4endl;
  }
#endif
  yEnd=  yCurrent;  
  return stepTrial; 
}


// ...........................................................................

G4double G4ChordFinder::NewStep(G4double  stepTrialOld, 
                                G4double  dChordStep, // Curr. dchord achieved
                                G4double& stepEstimate_Unconstrained )  
{
  // Is called to estimate the next step size, even for successful steps,
  // in order to predict an accurate 'chord-sensitive' first step
  // which is likely to assist in more performant 'stepping'.

  G4double stepTrial;

#if 1

  if (dChordStep > 0.0)
  {
    stepEstimate_Unconstrained =
                 stepTrialOld*std::sqrt( fDeltaChord / dChordStep );
    stepTrial =  fFractionNextEstimate * stepEstimate_Unconstrained;
  }
  else
  {
    // Should not update the Unconstrained Step estimate: incorrect!
    stepTrial =  stepTrialOld * 2.; 
  }

  if( stepTrial <= 0.001 * stepTrialOld)
  {
     if ( dChordStep > 1000.0 * fDeltaChord )
     {
        stepTrial= stepTrialOld * 0.03;   
     }
     else
     {
        if ( dChordStep > 100. * fDeltaChord )
        {
          stepTrial= stepTrialOld * 0.1;   
        }
        else   // Try halving the length until dChordStep OK
        {
          stepTrial= stepTrialOld * 0.5;   
        }
     }
  }
  else if (stepTrial > 1000.0 * stepTrialOld)
  {
     stepTrial= 1000.0 * stepTrialOld;
  }

  if( stepTrial == 0.0 )
  {
     stepTrial= 0.000001;
  }

#else

  if ( dChordStep > 1000. * fDeltaChord )
  {
        stepTrial= stepTrialOld * 0.03;   
  }
  else
  {
     if ( dChordStep > 100. * fDeltaChord )
     {
        stepTrial= stepTrialOld * 0.1;   
     }
     else  // Keep halving the length until dChordStep OK
     {
        stepTrial= stepTrialOld * 0.5;   
     }
  }

#endif 

  // A more sophisticated chord-finder could figure out a better
  // stepTrial, from dChordStep and the required d_geometry
  //   e.g.
  //      Calculate R, r_helix (eg at orig point)
  //      if( stepTrial < 2 pi  R )
  //          stepTrial = R arc_cos( 1 - fDeltaChord / r_helix )
  //      else    
  //          ??

  return stepTrial;
}


// ...........................................................................

G4FieldTrack
G4ChordFinder::ApproxCurvePointS( const G4FieldTrack&  CurveA_PointVelocity, 
                                  const G4FieldTrack&  CurveB_PointVelocity, 
                                  const G4FieldTrack&  ApproxCurveV,
                                  const G4ThreeVector& CurrentE_Point,
                                  const G4ThreeVector& CurrentF_Point,
                                  const G4ThreeVector& PointG,
                                       G4bool first, G4double eps_step)
{
  // ApproxCurvePointS is 2nd implementation of ApproxCurvePoint.
  // Use Brent Algorithm (or InvParabolic) when possible.
  // Given a starting curve point A (CurveA_PointVelocity), curve point B
  // (CurveB_PointVelocity), a point E which is (generally) not on the curve
  // and  a point F which is on the curve (first approximation), find new
  // point S on the curve closer to point E. 
  // While advancing towards S utilise 'eps_step' as a measure of the
  // relative accuracy of each Step.

  G4FieldTrack EndPoint(CurveA_PointVelocity);
  if(!first){EndPoint= ApproxCurveV;}

  G4ThreeVector Point_A,Point_B;
  Point_A=CurveA_PointVelocity.GetPosition();
  Point_B=CurveB_PointVelocity.GetPosition();

  G4double xa,xb,xc,ya,yb,yc;
 
  // InverseParabolic. AF Intersects (First Part of Curve) 

  if(first)
  {
    xa=0.;
    ya=(PointG-Point_A).mag();
    xb=(Point_A-CurrentF_Point).mag();
    yb=-(PointG-CurrentF_Point).mag();
    xc=(Point_A-Point_B).mag();
    yc=-(CurrentE_Point-Point_B).mag();
  }    
  else
  {
    xa=0.;
    ya=(Point_A-CurrentE_Point).mag();
    xb=(Point_A-CurrentF_Point).mag();
    yb=(PointG-CurrentF_Point).mag();
    xc=(Point_A-Point_B).mag();
    yc=-(Point_B-PointG).mag();
    if(xb==0.)
    {
      EndPoint=
      ApproxCurvePointV(CurveA_PointVelocity, CurveB_PointVelocity,
                        CurrentE_Point, eps_step);
      return EndPoint;
    }
  }

  const G4double tolerance= 1.e-12;
  if(std::abs(ya)<=tolerance||std::abs(yc)<=tolerance)
  {
    ; // What to do for the moment: return the same point as at start
      // then PropagatorInField will take care
  }
  else
  {
    G4double test_step = InvParabolic(xa,ya,xb,yb,xc,yc);
    G4double curve;
    if(first)
    {
      curve=std::abs(EndPoint.GetCurveLength()
                    -ApproxCurveV.GetCurveLength());
    }
    else
    {
      test_step=(test_step-xb);
      curve=std::abs(EndPoint.GetCurveLength()
                    -CurveB_PointVelocity.GetCurveLength());
      xb=(CurrentF_Point-Point_B).mag();
    }
      
    if(test_step<=0)    { test_step=0.1*xb; }
    if(test_step>=xb)   { test_step=0.5*xb; }
    if(test_step>=curve){ test_step=0.5*curve; } 

    if(curve*(1.+eps_step)<xb) // Similar to ReEstimate Step from
    {                          // G4VIntersectionLocator
      test_step=0.5*curve;
    }

    fIntgrDriver->AccurateAdvance(EndPoint,test_step, eps_step);
      
#ifdef G4DEBUG_FIELD
    // Printing Brent and Linear Approximation
    //
    G4cout << "G4ChordFinder::ApproxCurvePointS() - test-step ShF = "
           << test_step << "  EndPoint = " << EndPoint << G4endl;

    //  Test Track
    //
    G4FieldTrack TestTrack( CurveA_PointVelocity);
    TestTrack = ApproxCurvePointV( CurveA_PointVelocity, 
                                   CurveB_PointVelocity, 
                                   CurrentE_Point, eps_step );
    G4cout.precision(14);
    G4cout << "G4ChordFinder::BrentApprox = " << EndPoint  << G4endl;
    G4cout << "G4ChordFinder::LinearApprox= " << TestTrack << G4endl; 
#endif
  }
  return EndPoint;
}


// ...........................................................................

G4FieldTrack G4ChordFinder::
ApproxCurvePointV( const G4FieldTrack& CurveA_PointVelocity, 
                   const G4FieldTrack& CurveB_PointVelocity, 
                   const G4ThreeVector& CurrentE_Point,
                         G4double eps_step)
{
  // If r=|AE|/|AB|, and s=true path lenght (AB)
  // return the point that is r*s along the curve!
 
  G4FieldTrack   Current_PointVelocity = CurveA_PointVelocity; 

  G4ThreeVector  CurveA_Point= CurveA_PointVelocity.GetPosition();
  G4ThreeVector  CurveB_Point= CurveB_PointVelocity.GetPosition();

  G4ThreeVector  ChordAB_Vector= CurveB_Point   - CurveA_Point;
  G4ThreeVector  ChordAE_Vector= CurrentE_Point - CurveA_Point;

  G4double       ABdist= ChordAB_Vector.mag();
  G4double  curve_length;  //  A curve length  of AB
  G4double  AE_fraction; 
  
  curve_length= CurveB_PointVelocity.GetCurveLength()
              - CurveA_PointVelocity.GetCurveLength();  
 
  G4double  integrationInaccuracyLimit= std::max( perMillion, 0.5*eps_step ); 
  if( curve_length < ABdist * (1. - integrationInaccuracyLimit) )
  { 
#ifdef G4DEBUG_FIELD
    G4cerr << " Warning in G4ChordFinder::ApproxCurvePoint: "
           << G4endl
           << " The two points are further apart than the curve length "
           << G4endl
           << " Dist = "         << ABdist
           << " curve length = " << curve_length 
           << " relativeDiff = " << (curve_length-ABdist)/ABdist 
           << G4endl;
    if( curve_length < ABdist * (1. - 10*eps_step) )
    {
      std::ostringstream message;
      message << "Unphysical curve length." << G4endl
              << "The size of the above difference exceeds allowed limits."
              << G4endl
              << "Aborting.";
      G4Exception("G4ChordFinder::ApproxCurvePointV()", "GeomField0003",
                  FatalException, message);
    }
#endif
    // Take default corrective action: adjust the maximum curve length. 
    // NOTE: this case only happens for relatively straight paths.
    // curve_length = ABdist; 
  }

  G4double  new_st_length; 

  if ( ABdist > 0.0 )
  {
     AE_fraction = ChordAE_Vector.mag() / ABdist;
  }
  else
  {
     AE_fraction = 0.5;                         // Guess .. ?; 
#ifdef G4DEBUG_FIELD
     G4cout << "Warning in G4ChordFinder::ApproxCurvePointV():"
            << " A and B are the same point!" << G4endl
            << " Chord AB length = " << ChordAE_Vector.mag() << G4endl
            << G4endl;
#endif
  }
  
  if( (AE_fraction> 1.0 + perMillion) || (AE_fraction< 0.) )
  {
#ifdef G4DEBUG_FIELD
    G4cerr << " G4ChordFinder::ApproxCurvePointV() - Warning:"
           << " Anomalous condition:AE > AB or AE/AB <= 0 " << G4endl
           << "   AE_fraction = " <<  AE_fraction << G4endl
           << "   Chord AE length = " << ChordAE_Vector.mag() << G4endl
           << "   Chord AB length = " << ABdist << G4endl << G4endl;
    G4cerr << " OK if this condition occurs after a recalculation of 'B'"
           << G4endl << " Otherwise it is an error. " << G4endl ; 
#endif
     // This course can now result if B has been re-evaluated, 
     // without E being recomputed (1 July 99).
     // In this case this is not a "real error" - but it is undesired
     // and we cope with it by a default corrective action ...
     //
     AE_fraction = 0.5;                         // Default value
  }

  new_st_length= AE_fraction * curve_length; 

  if ( AE_fraction > 0.0 )
  { 
     fIntgrDriver->AccurateAdvance(Current_PointVelocity, 
                                   new_st_length, eps_step );
     //
     // In this case it does not matter if it cannot advance the full distance
  }

  // If there was a memory of the step_length actually required at the start 
  // of the integration Step, this could be re-used ...

  G4cout.precision(14);

  return Current_PointVelocity;
}


// ......................................................................

void
G4ChordFinder::PrintStatistics()
{
  // Print Statistics

  G4cout << "G4ChordFinder statistics report: " << G4endl;
  G4cout 
    << "  No trials: " << fTotalNoTrials_FNC
    << "  No Calls: "  << fNoCalls_FNC
    << "  Max-trial: " <<  fmaxTrials_FNC
    << G4endl; 
  G4cout 
    << "  Parameters: " 
    << "  fFirstFraction "  << fFirstFraction
    << "  fFractionLast "   << fFractionLast
    << "  fFractionNextEstimate " << fFractionNextEstimate
    << G4endl; 
}


// ...........................................................................

void G4ChordFinder::TestChordPrint( G4int    noTrials, 
                                    G4int    lastStepTrial, 
                                    G4double dChordStep, 
                                    G4double nextStepTrial )
{
     G4int oldprec= G4cout.precision(5);
     G4cout << " ChF/fnc: notrial " << std::setw( 3) << noTrials 
            << " this_step= "       << std::setw(10) << lastStepTrial;
     if( std::fabs( (dChordStep / fDeltaChord) - 1.0 ) < 0.001 )
     {
       G4cout.precision(8);
     }
     else
     {
       G4cout.precision(6);
     }
     G4cout << " dChordStep=  " << std::setw(12) << dChordStep;
     if( dChordStep > fDeltaChord ) { G4cout << " d+"; }
     else                           { G4cout << " d-"; }
     G4cout.precision(5);
     G4cout <<  " new_step= "       << std::setw(10)
            << fLastStepEstimate_Unconstrained
            << " new_step_constr= " << std::setw(10)
            << lastStepTrial << G4endl;
     G4cout << " nextStepTrial = " << std::setw(10) << nextStepTrial << G4endl;
     G4cout.precision(oldprec);
}
