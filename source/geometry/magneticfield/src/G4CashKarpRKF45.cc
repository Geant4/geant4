// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4CashKarpRKF45.cc,v 1.5 2000-11-01 15:15:52 gcosmo Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// The Cash-Karp Runge-Kutta-Fehlberg 4/5 method is an embedded fourth
//  order method (giving fifth-order accuracy) for the solution
//  of an ODE. Two different fourth order estimates are calculated;
//  their difference gives an error estimate.
// (We use it to integrate the equations of the motion of a particle 
//  in a magnetic field. )
//
//  Similar to Numerical Recipes, .... put REFerence here!
//

#include "G4CashKarpRKF45.hh"

/////////////////////////////////////////////////////////////////////
//
// Constructor

G4CashKarpRKF45::G4CashKarpRKF45(G4Mag_EqRhs *EqRhs, G4int numberOfVariables)
  : G4MagIntegratorStepper(EqRhs, numberOfVariables)
{
  fNumberOfVariables = numberOfVariables ;

  ak2 = new G4double[fNumberOfVariables] ;  
  ak3 = new G4double[fNumberOfVariables] ; 
  ak4 = new G4double[fNumberOfVariables] ; 
  ak5 = new G4double[fNumberOfVariables] ; 
  ak6 = new G4double[fNumberOfVariables] ; 
  yTemp = new G4double[fNumberOfVariables] ; 
  yIn = new G4double[fNumberOfVariables] ;
}

/////////////////////////////////////////////////////////////////////
//
// Destructor

G4CashKarpRKF45::~G4CashKarpRKF45()
{
  delete[] ak2;
  delete[] ak3;
  delete[] ak4;
  delete[] ak5;
  delete[] ak6;
  delete[] ak7;
  delete[] yTemp;
  delete[] yIn;
}

//////////////////////////////////////////////////////////////////////
//
// Given values for n = 6 variables yIn[0,...,n-1] 
// known  at x, use the fifth-order Cash-Karp Runge-
// Kutta-Fehlberg-4-5 method to advance the solution over an interval
// Step and return the incremented variables as yOut[0,...,n-1]. Also
// return an estimate of the local truncation error yErr[] using the
// embedded 4th-order method. The user supplies routine
// RightHandSide(y,dydx), which returns derivatives dydx for y .

void
G4CashKarpRKF45::Stepper(const G4double yInput[],
			 const G4double dydx[],
			       G4double Step,
			       G4double yOut[],
			       G4double yErr[])
{
  // const G4int nvar = 6 ;
 G4int i;

 const G4double a2 = 0.2 , a3 = 0.3 , a4 = 0.6 , a5 = 1.0 , a6 = 0.875 ,

                 b21 = 0.2 ,
                 b31 = 3.0/40.0 , b32 = 9.0/40.0 ,
                 b41 = 0.3 , b42 = -0.9 , b43 = 1.2 ,

                 b51 = -11.0/54.0 , b52 = 2.5 , b53 = -70.0/27.0 ,
                 b54 = 35.0/27.0 ,

                 b61 = 1631.0/55296.0 , b62 =   175.0/512.0 ,
                 b63 =  575.0/13824.0 , b64 = 44275.0/110592.0 ,
                 b65 =  253.0/4096.0 ,

                 c1 = 37.0/378.0 , c3 = 250.0/621.0 , c4 = 125.0/594.0 ,
                 c6 = 512.0/1771.0 ,
                                          dc5 = -277.0/14336.0 ;

 const G4double dc1 = c1 - 2825.0/27648.0 ,  dc3 = c3 - 18575.0/48384.0 ,
		dc4 = c4 - 13525.0/55296.0 , dc6 = c6 - 0.25 ;


   //  Saving yInput because yInput and yOut can be aliases for same array

   for(i=0;i<fNumberOfVariables;i++) 
   {
     yIn[i]=yInput[i];
   }
 // RightHandSide(yIn, dydx) ;              // 1st Step

 for(i=0;i<fNumberOfVariables;i++) 
 {
   yTemp[i] = yIn[i] + b21*Step*dydx[i] ;
 }
 RightHandSide(yTemp, ak2) ;              // 2nd Step

 for(i=0;i<fNumberOfVariables;i++)
 {
    yTemp[i] = yIn[i] + Step*(b31*dydx[i] + b32*ak2[i]) ;
 }
 RightHandSide(yTemp, ak3) ;              // 3rd Step

 for(i=0;i<fNumberOfVariables;i++)
 {
    yTemp[i] = yIn[i] + Step*(b41*dydx[i] + b42*ak2[i] + b43*ak3[i]) ;
 }
 RightHandSide(yTemp, ak4) ;              // 4th Step

 for(i=0;i<fNumberOfVariables;i++)
 {
    yTemp[i] = yIn[i] + Step*(b51*dydx[i] + b52*ak2[i] + b53*ak3[i] +
                      b54*ak4[i]) ;
 }
 RightHandSide(yTemp, ak5) ;              // 5th Step

 for(i=0;i<fNumberOfVariables;i++)
 {
    yTemp[i] = yIn[i] + Step*(b61*dydx[i] + b62*ak2[i] + b63*ak3[i] +
                      b64*ak4[i] + b65*ak5[i]) ;
 }
 RightHandSide(yTemp, ak6) ;              // 6th Step

 for(i=0;i<fNumberOfVariables;i++)
 {
    // Accumulate increments with proper weights

    yOut[i] = yIn[i] + Step*(c1*dydx[i] + c3*ak3[i] + c4*ak4[i] + c6*ak6[i]) ;

    // Estimate error as difference between 4th and
    // 5th order methods

    yErr[i] = Step*(dc1*dydx[i] + dc3*ak3[i] + dc4*ak4[i] +
              dc5*ak5[i] + dc6*ak6[i]) ;
 }
 // NormaliseTangentVector( yOut ); // Not wanted

 return ;

}   // end of Stepper .......................................................

void
G4CashKarpRKF45::StepWithEst(const G4double yInput[],
			     const G4double dydx[],
			           G4double Step,
			           G4double yOut[],
                                   G4double& alpha2,
			           G4double& beta2,
			     const G4double B1[],
			           G4double B2[]    )    
{

 G4Exception("G4CashKarpRKF45::StepWithEst ERROR: This Method is no longer used.");

#if 0
  // const G4int nvar = 6 ;
 G4int i;

   //  Saving yInput because yInput and yOut can be aliases for same array

   for(i=0;i<fNumberOfVariables;i++) 
   {
     yIn[i]=yInput[i];
   }

 const G4double a2 = 0.2 , a3 = 0.3 , a4 = 0.6 , a5 = 1.0 , a6 = 0.875 ,

                 b21 = 0.2 ,
                 b31 = 3.0/40.0 , b32 = 9.0/40.0 ,
                 b41 = 0.3 , b42 = -0.9 , b43 = 1.2 ,

                 b51 = -11.0/54.0 , b52 = 2.5 , b53 = -70.0/27.0 ,
                 b54 = 35.0/27.0 ,

                 b61 = 1631.0/55296.0 , b62 =   175.0/512.0 ,
                 b63 =  575.0/13824.0 , b64 = 44275.0/110592.0 ,
                 b65 =  253.0/4096.0 ,

                 c1 = 37.0/378.0 , c3 = 250.0/621.0 , c4 = 125.0/594.0 ,
                 c6 = 512.0/1771.0 ,
                                          dc5 = -277.0/14336.0 ;

 const G4double dc1 = c1 - 2825.0/27648.0 ,  dc3 = c3 - 18575.0/48384.0 ,
		dc4 = c4 - 13525.0/55296.0 , dc6 = c6 - 0.25 ;

 alpha2 = 0 ;
 beta2 = 0 ;

 // RightHandSide(yIn, dydx) ;              // 1st Step

 for(i=0;i<3;i++) 
 {
    alpha2 += B1[i]*B1[i] ;
    beta2 += dydx[i+3]*dydx[i+3] ;
 }

 for(i=0;i<fNumberOfVariables;i++)
 { 
   yTemp[i] = yIn[i] + b21*Step*dydx[i] ;
 }

 GetEquationOfMotion()->EvaluateRhsReturnB(yTemp,ak2,B2) ;  //  Calculates yderive & 
                                                    // returns B too!
 for(i=0;i<3;i++) 
 {
    alpha2 += B2[i]*B2[i] ;
    beta2 += ak2[i+3]*ak2[i+3] ;
 }

 for(i=0;i<fNumberOfVariables;i++)
 {
    yTemp[i] = yIn[i] + Step*(b31*dydx[i] + b32*ak2[i]) ;
 }
 GetEquationOfMotion()->EvaluateRhsReturnB(yTemp,ak3,B2) ;   
 for(i=0;i<3;i++)
 {
    alpha2 += B2[i]*B2[i] ;
    beta2 += ak3[i+3]*ak3[i+3] ;
 }

 for(i=0;i<fNumberOfVariables;i++)
 {
    yTemp[i] = yIn[i] + Step*(b41*dydx[i] + b42*ak2[i] + b43*ak3[i]) ;
 }
 GetEquationOfMotion()->EvaluateRhsReturnB(yTemp,ak4,B2) ;  
 
 for(i=0;i<3;i++)
 {
    alpha2 += B2[i]*B2[i] ;
    beta2 += ak4[i+3]*ak4[i+3] ;
 }

 for(i=0;i<fNumberOfVariables;i++)
 {
    yTemp[i] = yIn[i] + Step*(b51*dydx[i] + b52*ak2[i] + b53*ak3[i] +
                      b54*ak4[i]) ;
 }
 GetEquationOfMotion()->EvaluateRhsReturnB(yTemp,ak5,B2) ;
   
 for(i=0;i<3;i++)
 {
    alpha2 += B2[i]*B2[i] ;
    beta2 += ak5[i+3]*ak5[i+3] ;
 }

 for(i=0;i<fNumberOfVariables;i++)
 {
    yTemp[i] = yIn[i] + Step*(b61*dydx[i] + b62*ak2[i] + b63*ak3[i] +
                      b64*ak4[i] + b65*ak5[i]) ;
 }
 GetEquationOfMotion()->EvaluateRhsReturnB(yTemp,ak6,B2) ;  
 
 for(i=0;i<3;i++) 
 {
    alpha2 += B2[i]*B2[i] ;
    beta2 += ak6[i+3]*ak6[i+3] ;
 }

 for(i=0;i<fNumberOfVariables;i++)
 {
    // Accumulate increments with proper weights

    yOut[i] = yIn[i] + Step*(c1*dydx[i] + c3*ak3[i] + c4*ak4[i] + c6*ak6[i]) ;
 }

 alpha2 *= sqr(GetEquationOfMotion()->FCof()*Step)/6.0  ;
 beta2 *= (Step*Step)/6.0 ; 
 // NormaliseTangentVector( yOut );
#endif

 return ;

}   // end of StepWithEst ....................................................

