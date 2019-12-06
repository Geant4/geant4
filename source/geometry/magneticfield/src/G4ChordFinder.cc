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
// G4ChordFinder implementation
//
// Author: J.Apostolakis - Design and implementation - 25.02.1997
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
#include "G4NystromRK4.hh"

// New FSAL type driver / steppers -----
#include "G4IntegrationDriver.hh"
#include "G4InterpolationDriver.hh"
// #include "G4FSALBogackiShampine45.hh"
// #include "G4FSALDormandPrince745.hh"
#include "G4HelixHeum.hh"
#include "G4BFieldIntegrationDriver.hh"

#include <cassert>

// ..........................................................................

G4ChordFinder::G4ChordFinder(G4VIntegrationDriver* pIntegrationDriver)
  : fDefaultDeltaChord(0.25 * mm), fIntgrDriver(pIntegrationDriver)
{
  // Simple constructor -- it does not create equation

  fDeltaChord = fDefaultDeltaChord;       // Parameters
}


// ..........................................................................

G4ChordFinder::G4ChordFinder( G4MagneticField*        theMagField,
                              G4double                stepMinimum, 
                              G4MagIntegratorStepper* pItsStepper,
                              G4bool                  useFSALstepper )
  : fDefaultDeltaChord(0.25 * mm)
{
  // Construct the Chord Finder
  // by creating in inverse order the Driver, the Stepper and EqRhs ...

  fDeltaChord = fDefaultDeltaChord;       // Parameters

  using NewFsalStepperType = G4RK547FEq1; // or 2 or 3
  const char* NewFSALStepperName =
      "G4RK574FEq1> FSAL 4th/5th order 7-stage 'Equilibrium-type' #1.";
  using RegularStepperType =
         G4DormandPrince745; // 5th order embedded method. High efficiency.
         // G4ClassicalRK4;        // The old default
         // G4CashKarpRKF45;       // First embedded method in G4
         // G4BogackiShampine45;   // High efficiency 5th order embedded method
         // G4NystromRK4;          // Nystrom stepper 4th order 
         // G4RK547FEq1;  // or 2 or 3 
  const char* RegularStepperName = 
      "G4DormandPrince745 (aka DOPRI5): 5th/4th Order 7-stage embedded stepper";
      // "BogackiShampine 45 (Embedded 5th/4th Order, 7-stage)";
      // "Nystrom stepper 4th order";

  // Configurable
  G4bool forceFSALstepper = false; //  Choice - true to enable !!
  G4bool recallFSALflag  = useFSALstepper;
  useFSALstepper   = forceFSALstepper || useFSALstepper;

#ifdef G4DEBUG_FIELD
     G4cout << "G4ChordFinder 2nd Constructor called. " << G4endl;
     G4cout << " Parameters: " << G4endl;
     G4cout << "    useFSAL stepper= " << useFSALstepper
            << " (request = " << recallFSALflag 
            << " force FSAL = " << forceFSALstepper << " )" << G4endl;
#endif

  // useHigherStepper = forceHigherEffiencyStepper || useHigherStepper;
  
  G4Mag_EqRhs* pEquation = new G4Mag_UsualEqRhs(theMagField);
  fEquation = pEquation;                            

  // G4MagIntegratorStepper* regularStepper = nullptr;
  // G4VFSALIntegrationStepper* fsalStepper = nullptr; // for FSAL steppers only
  // G4MagIntegratorStepper* oldFSALStepper = nullptr;

  G4bool errorInStepperCreation = false;

  std::ostringstream message;  // In case of failure, load with description !

  if( pItsStepper != nullptr )
  {
     // Type is not known - so must use old class
     fIntgrDriver = new G4IntegrationDriver<G4MagIntegratorStepper>(
        stepMinimum, pItsStepper, pItsStepper->GetNumberOfVariables());
  }
  else if ( !useFSALstepper )
  {
     // RegularStepperType* regularStepper = nullptr; // To check the exception
     auto regularStepper = new RegularStepperType(pEquation);
     //                    *** ******************
     //
     // Alternative - for G4NystromRK4:
     // = new G4NystromRK4(pEquation, 0.1*mm );
     fRegularStepperOwned = regularStepper;

     if( regularStepper == nullptr )
     {
        message << "Stepper instantiation FAILED." << G4endl;        
        message << "G4ChordFinder: Attempted to instantiate "
                << RegularStepperName << " type stepper " << G4endl;
        G4Exception("G4ChordFinder::G4ChordFinder()",
                    "GeomField1001", JustWarning, message);
        errorInStepperCreation = true;
     }
     else
     {
        using SmallStepDriver = G4InterpolationDriver<G4DormandPrince745>;
        using LargeStepDriver = G4IntegrationDriver<G4HelixHeum>;

        fLongStepper = std::unique_ptr<G4HelixHeum>(new G4HelixHeum(pEquation));
        
        fIntgrDriver = new G4BFieldIntegrationDriver(
          std::unique_ptr<SmallStepDriver>(new SmallStepDriver(stepMinimum,
              regularStepper, regularStepper->GetNumberOfVariables())),
          std::unique_ptr<LargeStepDriver>(new LargeStepDriver(stepMinimum,
              fLongStepper.get(), regularStepper->GetNumberOfVariables())) );
        
        if( fIntgrDriver == nullptr)
        {        
           message << "Using G4BFieldIntegrationDriver with "
                   << RegularStepperName << " type stepper " << G4endl;
           message << "Driver instantiation FAILED." << G4endl;
           G4Exception("G4ChordFinder::G4ChordFinder()",
                       "GeomField1001", JustWarning, message);
        }
     }
  }
  else
  {
     auto fsalStepper=  new NewFsalStepperType(pEquation);
     //                 *** ******************
     fNewFSALStepperOwned = fsalStepper;

     if( fsalStepper == nullptr )
     {
        message << "Stepper instantiation FAILED." << G4endl;        
        message << "Attempted to instantiate "
                << NewFSALStepperName << " type stepper " << G4endl;
        G4Exception("G4ChordFinder::G4ChordFinder()",
                    "GeomField1001", JustWarning, message);
        errorInStepperCreation = true;
     }
     else
     {
        fIntgrDriver = new
           G4FSALIntegrationDriver<NewFsalStepperType>(stepMinimum, fsalStepper,
                                          fsalStepper->GetNumberOfVariables() );
           //  ====  Create the driver which knows the class type
        
        if( fIntgrDriver == nullptr )
        {
           message << "Using G4FSALIntegrationDriver with stepper type: "
                   << NewFSALStepperName << G4endl;
           message << "Integration Driver instantiation FAILED." << G4endl;
           G4Exception("G4ChordFinder::G4ChordFinder()",
                       "GeomField1001", JustWarning, message);
        }
     }
  }

  // -- Main work is now done
  
  //    Now check that no error occured, and report it if one did.
  
  // To test failure to create driver
  // delete fIntgrDriver;
  // fIntgrDriver = nullptr;

  // Detect and report Error conditions
  //
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
        errmsg  << "ERROR> Failure to create Integration-Driver object."
                << G4endl
                << "       -------------------------------------------"
                << G4endl;
     }
     const std::string BoolName[2]= { "False", "True" }; 
     errmsg << "  Configuration:  (constructor arguments) " << G4endl        
            << "    provided Stepper = " << pItsStepper << G4endl
            << "    use FSAL stepper = " << BoolName[useFSALstepper]
            << " (request = " << BoolName[recallFSALflag]
            << " force FSAL = " << BoolName[forceFSALstepper] << " )"
            << G4endl;
     errmsg << message.str(); 
     errmsg << "Aborting.";
     G4Exception("G4ChordFinder::G4ChordFinder() - constructor 2",
                 "GeomField0003", FatalException, errmsg);     
  }

  assert(    ( pItsStepper != nullptr ) 
          || ( fRegularStepperOwned != nullptr )
          || ( fNewFSALStepperOwned != nullptr )
     );
  assert( fIntgrDriver != nullptr );
}


// ......................................................................

G4ChordFinder::~G4ChordFinder()
{
  delete fEquation;
  delete fRegularStepperOwned;
  delete fNewFSALStepperOwned;
  delete fCachedField;
  delete fIntgrDriver;
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
  if(!first) { EndPoint = ApproxCurveV; }

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
      EndPoint = ApproxCurvePointV(CurveA_PointVelocity, CurveB_PointVelocity,
                                   CurrentE_Point, eps_step);
      return EndPoint;
    }
  }

  const G4double tolerance = 1.e-12;
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
      test_step = test_step - xb;
      curve=std::abs(EndPoint.GetCurveLength()
                    -CurveB_PointVelocity.GetCurveLength());
      xb = (CurrentF_Point-Point_B).mag();
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
 
  G4double integrationInaccuracyLimit= std::max( perMillion, 0.5*eps_step ); 
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

  G4double new_st_length; 

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

  new_st_length = AE_fraction * curve_length; 

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

// ...........................................................................
