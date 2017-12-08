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
// $Id: G4ChordFinder.cc 107508 2017-11-20 08:23:14Z gcosmo $
//
//
// 25.02.97 - John Apostolakis - Design and implementation 
// -------------------------------------------------------------------

#include <iomanip>

#include "G4ChordFinder.hh"
#include "G4SystemOfUnits.hh"
#include "G4MagneticField.hh"
#include "G4Mag_UsualEqRhs.hh"
#include "G4MagIntegratorDriver.hh"
// #include "G4ClassicalRK4.hh"
// #include "G4CashKarpRKF45.hh"
// #include "G4BogackiShampine23.hh"
// #include "G4BogackiShampine45.hh"
#include "G4DormandPrince745.hh"

// New FSAL type driver / steppers -----
#include "G4FSALIntegrationDriver.hh"
#include "G4VFSALIntegrationStepper.hh"
#include "G4RK547FEq1.hh"
// #include "G4RK547FEq2.hh"
// #include "G4RK547FEq3.hh"
// #include "G4NystromRK4.hh"

// New FSAL type driver / steppers -----
#include "G4IntegrationDriver.hh"
#include "G4FSALBogackiShampine45.hh"
// #include "G4FSALDormandPrince745.hh"


// ..........................................................................

G4ChordFinder::G4ChordFinder(G4VIntegrationDriver* pIntegrationDriver)
  : fDefaultDeltaChord( 0.25 * mm ),      // Parameters
    fDeltaChord( fDefaultDeltaChord ),    //   Internal parameters
    fFirstFraction(0.999), fFractionLast(1.00),  fFractionNextEstimate(0.98), 
    fMultipleRadius(15.0), 
    fStatsVerbose(0),
    fRegularStepperOwned(nullptr),                    // Dependent objects 
    fEquation(0),      
    fTotalNoTrials_FNC(0), fNoCalls_FNC(0), fmaxTrials_FNC(0)
{
  // Simple constructor -- it does not create equation
  fIntgrDriver= pIntegrationDriver;

  fLastStepEstimate_Unconstrained = DBL_MAX;          // Should move q, p to

  SetFractions_Last_Next( fFractionLast, fFractionNextEstimate);  
    // check the values and set the other parameters

  // G4cout << "G4ChordFinder 1st Constructor called - (driver given). " << G4endl;
}


// ..........................................................................

G4ChordFinder::G4ChordFinder( G4MagneticField*        theMagField,
                              G4double                stepMinimum, 
                              G4MagIntegratorStepper* pItsStepper,     // nullptr is default
                              // G4bool      useHigherEfficiencyStepper,  // false by default
                              G4bool                  useFSALstepper ) // false by default
  : fDefaultDeltaChord( 0.25 * mm ),     // Constants 
    fDeltaChord( fDefaultDeltaChord ),   // Parameters
    fFirstFraction(0.999), fFractionLast(1.00),  fFractionNextEstimate(0.98), 
    fMultipleRadius(15.0), 
    fStatsVerbose(0),
    // fRegularStepperOwned(nullptr),                    // Dependent objects     
    fEquation(0), 
    fTotalNoTrials_FNC(0), fNoCalls_FNC(0), fmaxTrials_FNC(0)  // State - stats
{
  //  Construct the Chord Finder
  //  by creating in inverse order the  Driver, the Stepper and EqRhs ...

  using NewFsalStepperType = G4RK547FEq1; // or 2 or 3
  const char* NewFSALStepperName = "G4RK574FEq1> FSAL 4th/5th order 7-stage 'Equilibrium-type' #1.";
//  using OldFsalStepperType = G4FSALBogackiShampine45;
//  const char* OldFSALStepperName = "FSAL BogackiShampine 45 (Embedded 5th/4th Order, 7-stage)";
          // = G4FSALDormandPrince745; // = "FSAL Dormand Prince 745 stepper";   
  using RegularStepperType =
         G4DormandPrince745; // Famous DOPRI5 (MatLab) 5th order embedded method. High efficiency.
         // G4ClassicalRK4;        // The old default
         // G4CashKarpRKF45;       // First embedded method in G4
         // G4BogackiShampine45;   // High efficiency 5th order embedded method
         // G4NystromRK4(pEquation, 0.1*millimeter ); // *clhep::millimeter );
         // G4RK547FEq1;  // or 2 or 3 
  const char* RegularStepperName = "G4DormandPrince745 (aka DOPRI5): 5th/4th Order 7-stage embedded stepper";
      // "BogackiShampine 45 (Embedded 5th/4th Order, 7-stage)";

  // Configurable
  G4bool forceFSALstepper= false; //  Choice - true to enable !!
  // G4bool useNewFSALtype= true;
  // G4bool forceHigherEffiencyStepper = false;
  G4bool report = false;  // Report type of stepper used

  bool recallFSALflag  = useFSALstepper;
  useFSALstepper   = forceFSALstepper || useFSALstepper;

  if( report ) {
     G4cout << "G4ChordFinder 2nd Constructor called. " << G4endl;
     G4cout << " Parameters: " << G4endl;
     G4cout << "    useFSAL stepper= " << useFSALstepper
            << " (request = " << recallFSALflag 
            << " force FSAL = " << forceFSALstepper << " )" << G4endl;
  }

  // useHigherStepper = forceHigherEffiencyStepper || useHigherStepper;
  
  G4Mag_EqRhs *pEquation = new G4Mag_UsualEqRhs(theMagField);
  fEquation = pEquation;                            
  fLastStepEstimate_Unconstrained = DBL_MAX;          // Should move q, p to
                                                     //    G4FieldTrack ??
  SetFractions_Last_Next( fFractionLast, fFractionNextEstimate);  
    // check the values and set the other parameters

  // G4MagIntegratorStepper*    regularStepper = nullptr;
  // G4VFSALIntegrationStepper*    fsalSepper  = nullptr; // for new-type FSAL steppers only
  // NewFsalStepperType*           fsalStepper =nullptr;
  // G4MagIntegratorStepper*    oldFSALStepper =nullptr;

  G4bool errorInStepperCreation = false;

  std::ostringstream message;  // In case of failure, load with description !
  message << "G4ChordFinder 2nd Constructor called. " << G4endl;

  if( pItsStepper != nullptr )
  {
     // Type is not known - so must use old class
     fIntgrDriver = new G4MagInt_Driver(stepMinimum, pItsStepper, 
                                        pItsStepper->GetNumberOfVariables() );
  }
  else if ( !useFSALstepper )
  {
     // RegularStepperType* regularStepper =nullptr;  // To check the exception
     auto regularStepper =new RegularStepperType(pEquation);
     //                   *** ******************
     //
     // Alternative - for G4NystromRK4:
     // = new G4NystromRK4(pEquation, 0.1*millimeter ); // *clhep::millimeter );
     fRegularStepperOwned = regularStepper;

     if( regularStepper == nullptr )
     {
        message << "  ERROR> 'Regular' RK Stepper instantiation FAILED." << G4endl;        
        message << "G4ChordFinder: Attempted to instantiate "
                << RegularStepperName << " type stepper " << G4endl;
        errorInStepperCreation = true;
     }
     else
     {
        fIntgrDriver =
           new G4MagInt_Driver(stepMinimum,
                               regularStepper,
                               regularStepper->GetNumberOfVariables() );
           //  ====  Create the old type of driver

           // Alternative: 
           // new G4IntegrationDriver<RegularStepperType>(stepMinimum,
           //  ====  Create the driver which knows the class type
        
        if( (fIntgrDriver==nullptr) || report ) {        
           message << "G4ChordFinder: Using G4IntegrationDriver with "
                   << RegularStepperName << " type stepper " << G4endl;
        }
        if(fIntgrDriver==nullptr) {
           message << "  ERROR> 'Regular' RK Driver instantiation FAILED." << G4endl;
        }
     }
  }
  else
  {
     auto fsalStepper=  new NewFsalStepperType(pEquation);
     //                     ******************
     fNewFSALStepperOwned = fsalStepper;
     // delete fsalStepper;
     // /*NewFsalStepperType* */ fsalStepper =nullptr;  // To check the exception

     if( fsalStepper == nullptr )
     {
        message << "  ERROR> 'FSAL' RK Stepper instantiation FAILED." << G4endl;        
        message << "G4ChordFinder: Attempted to instantiate "
                << NewFSALStepperName << " type stepper " << G4endl;
        errorInStepperCreation = true;
     }
     else
     {
        fIntgrDriver = new
           G4FSALIntegrationDriver<NewFsalStepperType>(stepMinimum,
                                                       fsalStepper,
                                                       fsalStepper->GetNumberOfVariables() );
           //  ====  Create the driver which knows the class type
        
        if( (fIntgrDriver==nullptr) || report ) {
           message << "G4ChordFinder: Using G4FSALIntegrationDriver with stepper type: " << G4endl
                   << NewFSALStepperName << " (new-FSAL type stepper.) " << G4endl;
        }
        if(fIntgrDriver==nullptr) {
           message << "  ERROR> FSAL Integration Driver instantiation FAILED." << G4endl;
        }
     }
  }

  // -- Main work is now done
  
  //    Now check that no error occured, and report it if one did.
  
  // To test failure to create driver
  // delete fIntgrDriver;
  // fIntgrDriver= nullptr;

  // Detect and report Error conditions
  if( errorInStepperCreation || (fIntgrDriver == nullptr ))
  {
     std::ostringstream errmsg;
     
     if( errorInStepperCreation )
     {
        errmsg  << "ERROR> Failure to create Stepper object." << G4endl
                << "       --------------------------------" << G4endl;
     }
     if (fIntgrDriver == nullptr )
     {
        errmsg  << "ERROR> Failure to create Integration-Driver object." << G4endl
                << "       -------------------------------------------" << G4endl;
     }
     const std::string BoolName[2]= { "False", "True" }; 
     errmsg << "  Configuration:  (constructor arguments) " << G4endl        
            << "    provided Stepper = " << pItsStepper << G4endl
            << "    use FSAL stepper = " << BoolName[useFSALstepper]
                               // ( useFSALstepper    ? "True" : "False" )
            << " (request = " << BoolName[recallFSALflag]
            << " force FSAL = " << BoolName[forceFSALstepper] << " )"  << G4endl;
//          << "     use new FSAL stp = " << ( useNewFSALstepper ? "True" : "False" ) << G4endl;
     errmsg << message.str(); 
     errmsg << "Aborting.";
     G4Exception("G4ChordFinder::G4ChordFinder() - constructor 2", "GeomField0003",
                  FatalException, errmsg);     
  }
  else if ( report )
  {
     G4cout << message.str();
  }

  assert(    ( pItsStepper    != nullptr ) 
          || ( fRegularStepperOwned != nullptr )
          || ( fNewFSALStepperOwned  != nullptr )
 //       || ( fOldFSALStepperOwned  != nullptr )
     );
  assert( fIntgrDriver != nullptr );

}


// ......................................................................

G4ChordFinder::~G4ChordFinder()
{
  delete   fEquation; // fIntgrDriver->pIntStepper->theEquation_Rhs;
  delete   fRegularStepperOwned;
  delete   fNewFSALStepperOwned;
  // delete   fOldFSALStepperOwned;       
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
                                    const G4ThreeVector& latestSafetyOrigin,
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

  // G4cout << ">G4ChordFinder::FindNextChord called." << G4endl;
   
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
  const unsigned int  maxTrials= 75; // Avoid endless loop for bad convergence 

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
  while( (! validEndPoint) && (noTrials < maxTrials) );
  // Loop checking, 07.10.2016, J. Apostolakis

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
