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
// If Stepping Angle ( h / R_curve) < pi/3  use Classical RK4Stepper
// Else use  HelixExplicitEuler Stepper
//
// History: 
// Derived from ExactHelicalStepper 18/05/07
//
// -------------------------------------------------------------------------

#include "G4HelixMixedStepper.hh"
#include "G4ClassicalRK4.hh"
#include "G4ThreeVector.hh"
G4HelixMixedStepper::G4HelixMixedStepper(G4Mag_EqRhs *EqRhs)
  : G4MagHelicalStepper(EqRhs)
    
{
  
  fRK4Stepper= new G4ClassicalRK4(EqRhs);
}

G4HelixMixedStepper::~G4HelixMixedStepper() {
     delete(fRK4Stepper);
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
   fLastStepSize=Step;

   if(Ang_curve<0.33*pi){
    
    fRK4Stepper->Stepper(yInput,dydx,Step,yOut,yErr);

   }
    else{
      const G4int nvar = 6 ;
      G4int i;
      G4double      yTemp[7], yIn[7] ;
      G4ThreeVector  Bfld_midpoint;
    //  Saving yInput because yInput and yOut can be aliases for same array
        for(i=0;i<nvar;i++) yIn[i]=yInput[i];

      G4double h = Step * 0.5;
 
     // Do two half steps
          AdvanceHelix(yIn,   Bfld,  h, yTemp);
          MagFieldEvaluate(yTemp, Bfld_midpoint) ;     
          AdvanceHelix(yTemp, Bfld_midpoint, h, yOut);
     // Do a full step    
          h = Step ;
          AdvanceHelix(yIn, Bfld, h, yTemp); 
     // Error estimation
          for(i=0;i<nvar;i++) {
          yErr[i] = yOut[i] - yTemp[i] ;
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
// ---------------------------------------------------------------------------

G4double G4HelixMixedStepper::DistChord()   const 
{
  // Implementation : must check whether h/R > 2 pi  !!
  //   If( h/R <  pi) use G4LineSection::DistLine
  //   Else           DistChord=R_helix
  //
  G4double distChord;
  G4double H_helix;
  H_helix=fLastStepSize; 
  G4double Ang_curve=GetAngCurve();
  
  if(Ang_curve<pi){
    
    distChord=0.5*H_helix*std::tan(0.25*Ang_curve);  

  }
  else{
    distChord=GetRadHelix();
  }
 
  return distChord;
  
}
