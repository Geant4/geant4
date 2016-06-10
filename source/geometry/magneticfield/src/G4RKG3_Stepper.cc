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
// $Id: G4RKG3_Stepper.cc 68055 2013-03-13 14:43:28Z gcosmo $
//
// -------------------------------------------------------------------

#include "G4RKG3_Stepper.hh"
#include "G4LineSection.hh"
#include "G4Mag_EqRhs.hh"

G4RKG3_Stepper::G4RKG3_Stepper(G4Mag_EqRhs *EqRhs)
  : G4MagIntegratorStepper(EqRhs,6), hStep(0.)
{
}

G4RKG3_Stepper::~G4RKG3_Stepper()
{
}

void G4RKG3_Stepper::Stepper(  const G4double yInput[8],
                               const G4double dydx[6],
                                     G4double Step,
                                     G4double yOut[8],
                                     G4double yErr[])
{
   G4double  B[3];
   G4int nvar = 6 ;
   G4int i;
   G4double  by15 = 1. / 15. ; // was  0.066666666 ;

   G4double yTemp[8], dydxTemp[6], yIn[8] ;
   //  Saving yInput because yInput and yOut can be aliases for same array
   for(i=0;i<nvar;i++) yIn[i]=yInput[i];
   yIn[6] = yInput[6];
   yIn[7] = yInput[7];
   G4double h = Step * 0.5; 
   hStep=Step;
   // Do two half steps

   StepNoErr(yIn, dydx,h, yTemp,B) ;
   
   //Store Bfld for DistChord Calculation
   for(i=0;i<3;i++)BfldIn[i]=B[i];

   //   RightHandSide(yTemp,dydxTemp) ;

   GetEquationOfMotion()->EvaluateRhsGivenB(yTemp,B,dydxTemp) ;  
   StepNoErr(yTemp,dydxTemp,h,yOut,B);      
        
   // Store midpoint, chord calculation
                                 
   fyMidPoint = G4ThreeVector( yTemp[0],  yTemp[1],  yTemp[2]); 

   // Do a full Step

   h *= 2 ;
   StepNoErr(yIn,dydx,h,yTemp,B); 
   for(i=0;i<nvar;i++)
   {
      yErr[i] = yOut[i] - yTemp[i] ;
      yOut[i] += yErr[i]*by15 ;          // Provides 5th order of accuracy
   }

   //Store values for DistChord method

   fyInitial = G4ThreeVector( yIn[0],   yIn[1],   yIn[2]);
   fpInitial = G4ThreeVector( yIn[3],   yIn[4],   yIn[5]);
   fyFinal   = G4ThreeVector( yOut[0],  yOut[1],  yOut[2]); 
  
   // NormaliseTangentVector( yOut );  // Deleted
}

// ---------------------------------------------------------------------------

// Integrator for RK from G3 with evaluation of error in solution and delta
// geometry based on naive similarity with the case of uniform magnetic field.
// B1[3] is input  and is the first magnetic field values
// B2[3] is output and is the final magnetic field values.

void G4RKG3_Stepper::StepWithEst( const G4double*,
                                  const G4double*,
                                        G4double,
                                        G4double*,
                                        G4double&,
                                        G4double&,
                                  const G4double*,
                                        G4double* )
   
{
  G4Exception("G4RKG3_Stepper::StepWithEst()", "GeomField0001",
              FatalException, "Method no longer used.");
}

// -----------------------------------------------------------------


// Integrator RK Stepper from G3 with only two field evaluation per Step. 
// It is used in propagation initial Step by small substeps after solution 
// error and delta geometry considerations. B[3] is magnetic field which 
// is passed from substep to substep.

void G4RKG3_Stepper::StepNoErr(const G4double tIn[8],
                               const G4double dydx[6],
                                     G4double Step,
                                     G4double tOut[8],
                                     G4double B[3]      )     // const
   
{ 
  
  //  Copy and edit the routine above, to delete alpha2, beta2, ...
   G4double K1[7],K2[7],K3[7],K4[7] ;
   G4double tTemp[8], yderiv[6] ;

  // Need Momentum value to give correct values to the coefficients in equation
  // Integration on unit velocity, but  tIn[3,4,5] is momentum 
   G4double mom,inverse_mom;
   G4int i ;
   const G4double c1=0.5,c2=0.125,c3=1./6.;
  
   // GetEquationOfMotion()->EvaluateRhsReturnB(tIn,dydx,B1) ;
   // Correction for momentum not a velocity
   // Need the protection !!! must be not zero 
     mom=std::sqrt(tIn[3]*tIn[3]+tIn[4]*tIn[4]+tIn[5]*tIn[5]); 
     inverse_mom=1./mom;    
   for(i=0;i<3;i++)
   {
      K1[i] = Step * dydx[i+3]*inverse_mom;
      tTemp[i] = tIn[i] + Step*(c1*tIn[i+3]*inverse_mom + c2*K1[i]) ;
      tTemp[i+3] = tIn[i+3] + c1*K1[i]*mom ;
     
   }
    
   GetEquationOfMotion()->EvaluateRhsReturnB(tTemp,yderiv,B) ;
    
      
   for(i=0;i<3;i++)
   {
      K2[i] = Step * yderiv[i+3]*inverse_mom;
      tTemp[i+3] = tIn[i+3] + c1*K2[i]*mom ;
   }
   
   //  Given B, calculate yderiv !
   GetEquationOfMotion()->EvaluateRhsGivenB(tTemp,B,yderiv) ;  
 
   for(i=0;i<3;i++)
   {
      K3[i] = Step * yderiv[i+3]*inverse_mom;
      tTemp[i] = tIn[i] + Step*(tIn[i+3]*inverse_mom + c1*K3[i]) ;
      tTemp[i+3] = tIn[i+3] + K3[i]*mom ;
   }
   

   //  Calculates y-deriv(atives) & returns B too!
   GetEquationOfMotion()->EvaluateRhsReturnB(tTemp,yderiv,B) ;  
    

   for(i=0;i<3;i++)        // Output trajectory vector
   {
      K4[i] = Step * yderiv[i+3]*inverse_mom;
      tOut[i] = tIn[i] + Step*(tIn[i+3]*inverse_mom+ (K1[i] + K2[i] + K3[i])*c3) ;
      tOut[i+3] = tIn[i+3] + mom*(K1[i] + 2*K2[i] + 2*K3[i] +K4[i])*c3 ;
   }
   tOut[6] = tIn[6];
   tOut[7] = tIn[7];
   // NormaliseTangentVector( tOut );
  

}


// ---------------------------------------------------------------------------
 
 G4double G4RKG3_Stepper::DistChord()   const 
 {
   // Soon: must check whether h/R > 2 pi  !!
   //  Method below is good only for < 2 pi
   G4double distChord,distLine;
   
   if (fyInitial != fyFinal) {
      distLine= G4LineSection::Distline(fyMidPoint,fyInitial,fyFinal );
  
        distChord = distLine;
   }else{
      distChord = (fyMidPoint-fyInitial).mag();
   }
 
  
   return distChord;
   
 }

