// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4ClassicalRK4.cc,v 1.2 1999-12-15 14:49:49 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
#include "G4ClassicalRK4.hh"
#include "G4ThreeVector.hh"


//////////////////////////////////////////////////////////////////
//
// Constructor sets the number of variables (default = 6)

G4ClassicalRK4:: G4ClassicalRK4(G4Mag_EqRhs *EqRhs, G4int numberOfVariables): 
G4MagErrorStepper(EqRhs, numberOfVariables),
  fNumberOfVariables(numberOfVariables)
{
   dydxm = new G4double[fNumberOfVariables];
   dydxt = new G4double[fNumberOfVariables]; 
   yt    = new G4double[fNumberOfVariables]; 
}

////////////////////////////////////////////////////////////////
//
// Destructor

G4ClassicalRK4:: ~G4ClassicalRK4()
{
  delete[] dydxm;
  delete[] dydxt;
  delete[] yt;
}

//////////////////////////////////////////////////////////////////////
//
// Given values for the variables y[0,..,n-1] and their derivatives
// dydx[0,...,n-1] known at x, use the classical 4th Runge-Kutta
// method to advance the solution over an interval h and return the
// incremented variables as yout[0,...,n-1], which not be a distinct
// array from y. The user supplies the routine RightHandSide(x,y,dydx),
// which returns derivatives dydx at x. The source is routine rk4 from
// NRC p. 712-713 .

void
G4ClassicalRK4::DumbStepper( const G4double  yIn[],
			     const G4double  dydx[],
			     const G4double  h,
			 	   G4double  yOut[])
{
  const G4int nvar = fNumberOfVariables;   //  GetNumberOfVariables(); 
  G4int i;
  G4double  hh = h*0.5 , h6 = h/6.0  ;


  for(i=0;i<nvar;i++)
  {
    yt[i] = yIn[i] + hh*dydx[i] ;             // 1st Step K1=h*dydx
  }
  RightHandSide(yt,dydxt) ;                   // 2nd Step K2=h*dydxt

  for(i=0;i<nvar;i++)
  { 
    yt[i] = yIn[i] + hh*dydxt[i] ;
  }
  RightHandSide(yt,dydxm) ;                   // 3rd Step K3=h*dydxm

  for(i=0;i<nvar;i++)
  {
    yt[i]   = yIn[i] + h*dydxm[i] ;
    dydxm[i] += dydxt[i] ;                    // now dydxm=(K2+K3)/h
  }
  RightHandSide(yt,dydxt) ;                   // 4th Step K4=h*dydxt
 
  for(i=0;i<nvar;i++)    // Final RK4 output
  {
    yOut[i] = yIn[i] + h6*(dydx[i]+dydxt[i]+2.0*dydxm[i]); //+K1/6+K4/6+(K2+K3)/3
  }
  // NormaliseTangentVector( yOut );
  
  return ;

}  // end of DumbStepper ....................................................

////////////////////////////////////////////////////////////////////
//
//   

void
G4ClassicalRK4::StepWithEst( const G4double  yIn[],
			     const G4double  dydx[],
			     const G4double  h,
			 	   G4double  yOut[],
                                   G4double& alpha2,
                                   G4double& beta2,
			     const G4double  B1[],
			           G4double  B2[]         ) 
{

 G4Exception(" G4ClassicalRK4::StepWithEst ERROR: this Method is no longer used.");

#if 0  //  const G4int nvar = 6 ; 
  const G4int nvar = GetNumberOfVariables(); 
  G4int i;
  G4double  hh = h*0.5 , h6 = h/6.0 ;
  G4double B[3] ;
  alpha2 = 0 ;
  beta2 = 0 ;

  for(i=0;i<nvar;i++)
  {
    yt[i] = yIn[i] + hh*dydx[i] ;             // 1st Step K1=h*dydx
  }
  GetEquationOfMotion()->EvaluateRhsReturnB(yt,dydxt,B) ;  //  Calculates yderive & 
                                                      // returns B too!
   //  RightHandSide(yt,dydxt) ;                   // 2nd Step K2=h*dydxt

  for(i=0;i<nvar;i++)
  { 
    yt[i] = yIn[i] + hh*dydxt[i] ;
  }
   GetEquationOfMotion()->EvaluateRhsReturnB(yt,dydxm,B2) ;  
  //  RightHandSide(yt,dydxm) ;                   // 3rd Step K3=h*dydxm

  for(i=0;i<3;i++)
  {
     beta2 +=  dydx[i+3] *dydx[i+3] 
	     + dydxt[i+3]*dydxt[i+3] 
             + dydxm[i+3]*dydxm[i+3] ;
     alpha2 += B1[i]*B1[i] + B[i]*B[i] + B2[i]*B2[i] ;
  }
  for(i=0;i<nvar;i++)
  {
    yt[i]   = yIn[i] + h*dydxm[i] ;
    dydxm[i] += dydxt[i] ;                    // now dydxm=(K2+K3)/h
  }
   GetEquationOfMotion()->EvaluateRhsReturnB(yt,dydxt,B2) ;  
  //  RightHandSide(yt,dydxt) ;                   // 4th Step K4=h*dydxt
 
  for(i=0;i<nvar;i++)    // Final RK4 output
  {
    yOut[i] = yIn[i] + h6*(dydx[i]+dydxt[i]+2.0*dydxm[i]); //+K1/6+K4/6+(K2+K3)/3
  }
  for(i=0;i<3;i++)
  { 
     beta2 += dydxt[i+3]*dydxt[i+3] ;
     alpha2 += B2[i]*B2[i] ;
  }
  beta2  *= 0.25*h*h;
  alpha2 *= sqr(GetEquationOfMotion()->FCof()*h)*0.25 ;
  // NormaliseTangentVector( yOut );

#endif

  return ;

}  // end of StepWithEst ......................................................
   


