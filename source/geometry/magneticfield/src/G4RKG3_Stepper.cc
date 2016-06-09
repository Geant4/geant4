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
// $Id: G4RKG3_Stepper.cc,v 1.14 2007/05/18 12:44:28 tnikitin Exp $
// GEANT4 tag $Name: geant4-09-00 $
//
// -------------------------------------------------------------------

#include "G4RKG3_Stepper.hh"
#include "G4LineSection.hh"
#include "G4Mag_EqRhs.hh"

G4RKG3_Stepper::G4RKG3_Stepper(G4Mag_EqRhs *EqRhs)
  : G4MagIntegratorStepper(EqRhs,6)
{
  
  fPtrMagEqOfMot=EqRhs;

}

G4RKG3_Stepper::~G4RKG3_Stepper()
{
}

void G4RKG3_Stepper::Stepper(  const G4double  yInput[7],
                               const G4double dydx[7],
                                     G4double Step,
                                     G4double yOut[7],
                                     G4double yErr[])
{
   G4double  B[3];
   G4int nvar = 6 ;
   G4int i;
   G4double  by15 = 1. / 15. ; // was  0.066666666 ;

   G4double yTemp[7], dydxTemp[6], yIn[7] ;
   //  Saving yInput because yInput and yOut can be aliases for same array
   for(i=0;i<nvar;i++) yIn[i]=yInput[i];

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
   StepNoErr(yIn,dydx,h,yTemp,B); // ,beTemp2) ;
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
  G4Exception("G4RKG3_Stepper::StepWithEst()", "ObsoleteMethod",
              FatalException, "Method no longer used.");
}

// -----------------------------------------------------------------


// Integrator RK Stepper from G3 with only two field evaluation per Step. 
// It is used in propagation initial Step by small substeps after solution 
// error and delta geometry considerations. B[3] is magnetic field which 
// is passed from substep to substep.

void G4RKG3_Stepper::StepNoErr(const G4double tIn[7],
                               const G4double dydx[7],
                                     G4double Step,
                                     G4double tOut[7],
                                     G4double B[3]      )     // const
   
{ 
  
  //  Copy and edit the routine above, to delete alpha2, beta2, ...
   G4double K1[7],K2[7],K3[7],K4[7] ;
   G4double tTemp[7], yderiv[6] ;
  // Need Momentum value to give correct values to the coefficients in equation
  // Integration on unit velocity, but  tIn[3,4,5] is momentum 
   G4double mom;
   G4int i ;
  
#ifdef END_CODE_G3STEPPER
   G4Exception(" G4RKG3_Stepper::StepNoErr(): method to be no longer used.");
#else
  
   // GetEquationOfMotion()->EvaluateRhsReturnB(tIn,dydx,B1) ;
   // Correction for momentum not a velocity
   // Need the protection !!! must be not zero 
     mom=std::sqrt(tIn[3]*tIn[3]+tIn[4]*tIn[4]+tIn[5]*tIn[5]); 
         
   for(i=0;i<3;i++)
   {
      K1[i] = Step * dydx[i+3]/mom;
      tTemp[i] = tIn[i] + Step*(0.5*tIn[i+3]/mom + 0.125*K1[i]) ;
      tTemp[i+3] = tIn[i+3] + 0.5*K1[i]*mom ;
     
   }
    
   GetEquationOfMotion()->EvaluateRhsReturnB(tTemp,yderiv,B) ;
    
      
   for(i=0;i<3;i++)
   {
      K2[i] = Step * yderiv[i+3]/mom;
      tTemp[i+3] = tIn[i+3] + 0.5*K2[i]*mom ;
   }
   
   //  Given B, calculate yderiv !
   GetEquationOfMotion()->EvaluateRhsGivenB(tTemp,B,yderiv) ;  
 
   for(i=0;i<3;i++)
   {
      K3[i] = Step * yderiv[i+3]/mom;
      tTemp[i] = tIn[i] + Step*(tIn[i+3]/mom + 0.5*K3[i]) ;
      tTemp[i+3] = tIn[i+3] + K3[i]*mom ;
   }
   

   //  Calculates y-deriv(atives) & returns B too!
   GetEquationOfMotion()->EvaluateRhsReturnB(tTemp,yderiv,B) ;  
    

   for(i=0;i<3;i++)        // Output trajectory vector
   {
      K4[i] = Step * yderiv[i+3]/mom;
      tOut[i] = tIn[i] + Step*(tIn[i+3]/mom + (K1[i] + K2[i] + K3[i])/6.0) ;
      tOut[i+3] = tIn[i+3] + mom*(K1[i] + 2*K2[i] + 2*K3[i] +K4[i])/6.0 ;
   }
  
   // NormaliseTangentVector( tOut );
  
#endif
}

// ---------------------------------------------------------------------------

G4double G4RKG3_Stepper::DistChord()   const 
{
  // Implementation : must check whether h/R >  pi  !!
  //   If( h/R <  pi) use  G4LineSection::DistLine
  //   Else           use  DistChord=R_helix
  G4double distChord,distLine;

  //Calculation of R_helix and R_curv
  G4double R_helix;
  G4double R_curv;
  G4double H_helix=hStep;
  G4double Bmag=BfldIn.mag();
    
  G4ThreeVector initVelocity= fpInitial;
  G4double      velocityVal = initVelocity.mag();
  G4ThreeVector initTangent = (1.0/velocityVal) * initVelocity;  // .unit();  
    
  const G4double fUnitConstant = 0.299792458 * (GeV/(tesla*m));
  G4double particleCharge = fPtrMagEqOfMot->FCof() / (eplus*c_light); 
  G4double fCoefficient = (fUnitConstant ) * particleCharge;
  
  // for too small field there is no curvature
  
  if( Bmag>1e-12 ) {
  
    // Bnorm = Bfld.unit();
    G4ThreeVector Bnorm = (1.0/Bmag)*BfldIn;

    // calculate the direction of the force
   
    G4ThreeVector B_x_P = Bnorm.cross(initTangent);

    // parallel and perp vectors

    G4ThreeVector B_x_P_x_B = B_x_P.cross(Bnorm); 
    G4double ptan=B_x_P_x_B.dot(initVelocity);
     
     R_helix  =std::abs((  1./fCoefficient)* ptan/Bmag);
     R_curv=std::abs((  1./fCoefficient)* velocityVal/Bmag);

     // G4cout<<"Bfld="<<BfldIn<<" Momentum="<<velocityVal<<" DirectionM="<<initTangent<<G4endl;
     // G4cout<<"R_helix="<<R_helix/mm<<" mm R_curv="<<R_curv/mm<<" BxP="<< B_x_P<<" BxPxP="<<B_x_P_x_B<<G4endl;   
  }
  else{
   
    R_helix=0.;
    R_curv=H_helix;
  }
   
  // DistChord Calculation 
  if(std::abs(H_helix/R_curv)<pi){
    if (fyInitial != fyFinal) {
      distLine= G4LineSection::Distline(fyMidPoint,fyInitial,fyFinal );
      distChord = distLine;
    }else{
      distChord = (fyMidPoint-fyInitial).mag();
    }
  }
  else{
    distChord=R_helix;
  }
  // G4cout<<"distChord="<<distChord<<" hstep="<<H_helix<<"  Helix/R ="<<std::abs(H_helix/R_curv)<<" R_helix="<<R_helix<<G4endl; 
  return distChord;
  
}
