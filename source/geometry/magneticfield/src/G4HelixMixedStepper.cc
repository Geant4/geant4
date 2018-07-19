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
// class G4HelixMixedStepper
//
// Class description:
//
// G4HelixMixedStepper split the Method used for Integration in two:
//
// If Stepping Angle ( h / R_curve) < pi/3  
//        use Stepper for small step(ClassicalRK4 by default)
// Else use  HelixExplicitEuler Stepper
//
// History: 
// Derived from ExactHelicalStepper 18/05/07
//
// -------------------------------------------------------------------------

#include "G4HelixMixedStepper.hh"
#include "G4PhysicalConstants.hh"
#include "G4ClassicalRK4.hh"
#include "G4CashKarpRKF45.hh"
#include "G4SimpleRunge.hh"
#include "G4HelixImplicitEuler.hh"
#include "G4HelixExplicitEuler.hh"
#include "G4HelixSimpleRunge.hh"
#include "G4ExactHelixStepper.hh"
#include "G4ExplicitEuler.hh"
#include "G4ImplicitEuler.hh"
#include "G4SimpleHeum.hh"
#include "G4RKG3_Stepper.hh"
#include "G4NystromRK4.hh"
// Additional potential stepper 
#include "G4DormandPrince745.hh"
#include "G4BogackiShampine23.hh"
#include "G4BogackiShampine45.hh"
#include "G4TsitourasRK45.hh"

#include "G4ThreeVector.hh"
#include "G4LineSection.hh"

G4HelixMixedStepper::
G4HelixMixedStepper(G4Mag_EqRhs *EqRhs, G4int stepperNumber,
                    G4double angleThreshold)
  : G4MagHelicalStepper(EqRhs), fNumCallsRK4(0), fNumCallsHelix(0)
{
   SetVerbose(1);
   if( angleThreshold < 0.0 ){
     fAngle_threshold= 0.33*pi;
   }else{
     fAngle_threshold= angleThreshold;
   }

   if(stepperNumber<0)
     stepperNumber=4;  // Default is RK4   (original)      
     // stepperNumber=745;   // Default is DormandPrince745 (ie DoPri5)
     // stepperNumber=8;  // Default is CashKarp

   fStepperNumber = stepperNumber; // Store the choice
   fRK4Stepper =  SetupStepper(EqRhs, fStepperNumber);
}

G4HelixMixedStepper::~G4HelixMixedStepper()
{  
  delete(fRK4Stepper);
  if (fVerbose>0){ PrintCalls();};
}

void G4HelixMixedStepper::Stepper(  const G4double  yInput[7],
                               const G4double dydx[7],
                                     G4double Step,
                                     G4double yOut[7],
                                     G4double yErr[])
{
  //Estimation of the Stepping Angle
  G4ThreeVector Bfld;
  MagFieldEvaluate(yInput, Bfld);
  
  G4double Bmag = Bfld.mag();
  const G4double *pIn = yInput+3;
  G4ThreeVector initVelocity= G4ThreeVector( pIn[0], pIn[1], pIn[2]);
  G4double      velocityVal = initVelocity.mag();
  G4double R_1;
  G4double Ang_curve;
  
  R_1=std::abs(GetInverseCurve(velocityVal,Bmag));
  Ang_curve=R_1*Step;
  SetAngCurve(Ang_curve);
  SetCurve(std::abs(1/R_1));
  
  if(Ang_curve< fAngle_threshold){
    fNumCallsRK4++;
    fRK4Stepper->Stepper(yInput,dydx,Step,yOut,yErr);
  }
  else
  {
    fNumCallsHelix++;
    const G4int    nvar    = 6 ;
    const G4int    nvarMax = 8 ;    
    G4int          i;
    G4double       yTemp[nvarMax], yIn[nvarMax], yTemp2[nvarMax];
    G4ThreeVector  Bfld_midpoint;

    //  Saving yInput because yInput and yOut can be aliases for same array
    for(i=0;i<nvar;i++) yIn[i]=yInput[i];
    
    G4double halfS = Step * 0.5;
    // 1. Do first half step and full step
    AdvanceHelix(yIn, Bfld, halfS, yTemp, yTemp2); // yTemp2 for s=2*h (halfS)
    //**********

    MagFieldEvaluate(yTemp, Bfld_midpoint) ;

    // 2. Do second half step - with revised field
    // NOTE: Could avoid this call if  'Bfld_midpoint == Bfld'
    //       or diff 'almost' zero
    AdvanceHelix(yTemp, Bfld_midpoint, halfS, yOut);
    // Not requesting y at s=2*h (halfS)
    //**********
    
    // 3. Estimate the integration error
    //    should be (nearly) zero if Bfield= constant
    for(i=0;i<nvar;i++) {
      yErr[i] = yOut[i] - yTemp2[i] ;
    }
  }
}

void
G4HelixMixedStepper::DumbStepper( const G4double  yIn[],
				   G4ThreeVector   Bfld,
				   G4double        h,
				   G4double        yOut[])
{
  AdvanceHelix(yIn, Bfld, h, yOut);
}

G4double G4HelixMixedStepper::DistChord()   const
{
  // Implementation : must check whether h/R > 2 pi  !!
  //   If( h/R <  pi) use G4LineSection::DistLine
  //   Else           DistChord=R_helix
  //
  G4double distChord;
  G4double Ang_curve=GetAngCurve();
  
  if(Ang_curve<=pi){
    distChord=GetRadHelix()*(1-std::cos(0.5*Ang_curve));
  }
  else
  {
    if(Ang_curve<twopi){
      distChord=GetRadHelix()*(1+std::cos(0.5*(twopi-Ang_curve)));
    }
    else{
      distChord=2.*GetRadHelix();
    }
  }
  
  return distChord;
}

// ---------------------------------------------------------------------------
void G4HelixMixedStepper::PrintCalls()
{
  G4cout << "In HelixMixedStepper::Number of calls to smallStepStepper = "
         << fNumCallsRK4
         << "  and Number of calls to Helix = " << fNumCallsHelix << G4endl;
}

G4MagIntegratorStepper*
G4HelixMixedStepper::SetupStepper(G4Mag_EqRhs* pE, G4int StepperNumber)
{
  G4MagIntegratorStepper* pStepper;
  if (fVerbose>0) G4cout << " G4HelixMixedStepper: ";
  switch ( StepperNumber )
    {
      // Robust, classic method
      case 4:
        pStepper = new G4ClassicalRK4( pE );
        if (fVerbose>0) G4cout << "G4ClassicalRK4";
        break;

      // Steppers with embedded estimation of error
      case 8:
        pStepper = new G4CashKarpRKF45( pE );
        if (fVerbose>0) G4cout << "G4CashKarpRKF45";
        break;
      case 13:
        pStepper = new G4NystromRK4( pE );
        if (fVerbose>0) G4cout << "G4NystromRK4";
        break;
        
      // Lowest order RK Stepper - experimental
      case 1:
        pStepper = new G4ImplicitEuler( pE );
        if (fVerbose>0) G4cout << "G4ImplicitEuler";
        break;

      // Lower order RK Steppers - ok overall, good for uneven fields        
      case 2:
        pStepper = new G4SimpleRunge( pE );
        if (fVerbose>0) G4cout << "G4SimpleRunge";
        break;
      case 3:
        pStepper = new G4SimpleHeum( pE );
        if (fVerbose>0) G4cout << "G4SimpleHeum";
        break;
      case 23:
        pStepper = new G4BogackiShampine23( pE );
        if (fVerbose>0) G4cout << "G4BogackiShampine23";
        break;

      // Higher order RK Steppers
      // for smoother fields and high accuracy requirements 
      case 45:
        pStepper = new G4BogackiShampine45( pE );
        if (fVerbose>0) G4cout << "G4BogackiShampine45";
        break;
      case 145:
        pStepper = new G4TsitourasRK45( pE );
        if (fVerbose>0) G4cout << "G4TsitourasRK45";
        break;
      case 745:
        pStepper = new G4DormandPrince745( pE );
        if (fVerbose>0) G4cout << "G4DormandPrince745";
        break;

      // Helical Steppers
      case 6:
        pStepper = new G4HelixImplicitEuler( pE );
        if (fVerbose>0) G4cout << "G4HelixImplicitEuler";
        break;
      case 7:
        pStepper = new G4HelixSimpleRunge( pE );
        if (fVerbose>0) G4cout << "G4HelixSimpleRunge";
        break;
      case 5:
        pStepper = new G4HelixExplicitEuler( pE );
        if (fVerbose>0) G4cout << "G4HelixExplicitEuler";
        break; //  Since Helix Explicit is used for long steps,
               // this is useful only to measure overhead.
      // Exact Helix - likely good only for cases of
      //            i) uniform field (potentially over small distances)
      //           ii) segmented uniform field (maybe)
      case 9:
        pStepper = new G4ExactHelixStepper( pE );
        if (fVerbose>0) G4cout << "G4ExactHelixStepper";
        break;
      case 10:
        pStepper = new G4RKG3_Stepper( pE );
        if (fVerbose>0) G4cout << "G4RKG3_Stepper";
        break;

      // Low Order Steppers - not good except for very weak fields
      case 11:
        pStepper = new G4ExplicitEuler( pE );
        if (fVerbose>0) G4cout << "G4ExplicitEuler";
        break;
      case 12:
        pStepper = new G4ImplicitEuler( pE );
        if (fVerbose>0) G4cout << "G4ImplicitEuler";
        break;

      case 0:
      case -1:
      default:
        pStepper = new G4ClassicalRK4( pE );
        if (fVerbose>0) G4cout << "G4ClassicalRK4 (Default)";
        break;
    }
  if(fVerbose>0)
    G4cout << " chosen as stepper for small steps in G4HelixMixedStepper."
           << G4endl;

  return pStepper;
}
