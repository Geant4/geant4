// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4GaussJacobiQ.cc,v 1.1 1999-01-07 16:08:55 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
#include "G4GaussJacobiQ.hh"


// -------------------------------------------------------------
//
// Constructor for Gauss-Jacobi integration method. 
//

G4GaussJacobiQ::G4GaussJacobiQ(       function pFunction,
			              G4double alpha,
                                      G4double beta, 
			              G4int nJacobi           ) 
   : G4VGaussianQuadrature(pFunction)

{
  const G4double tolerance = 1.0e-12 ;
  const G4double maxNumber = 12 ;
  G4int i, k, j ;
  G4double alphaBeta, alphaReduced, betaReduced, root1, root2, root3 ;
  G4double a, b, c, newton1, newton2, newton3, newton, temp, root, rootTemp ;

  fNumber   = nJacobi ;
  fAbscissa = new G4double[fNumber] ;
  fWeight   = new G4double[fNumber] ;

  for (i=1;i<=nJacobi;i++)
  {
     if (i == 1)
     {
	alphaReduced = alpha/nJacobi ;
	betaReduced = beta/nJacobi ;
	root1 = (1.0+alpha)*(2.78002/(4.0+nJacobi*nJacobi)+
	      0.767999*alphaReduced/nJacobi) ;
	root2 = 1.0+1.48*alphaReduced+0.96002*betaReduced +
	     0.451998*alphaReduced*alphaReduced+0.83001*alphaReduced*betaReduced ;
	root  = 1.0-root1/root2 ;
     } 
     else if (i == 2)
     {
	root1=(4.1002+alpha)/((1.0+alpha)*(1.0+0.155998*alpha)) ;
	root2=1.0+0.06*(nJacobi-8.0)*(1.0+0.12*alpha)/nJacobi ;
	root3=1.0+0.012002*beta*(1.0+0.24997*fabs(alpha))/nJacobi ;
	root -= (1.0-root)*root1*root2*root3 ;
     } 
     else if (i == 3) 
     {
	root1=(1.67001+0.27998*alpha)/(1.0+0.37002*alpha) ;
	root2=1.0+0.22*(nJacobi-8.0)/nJacobi ;
	root3=1.0+8.0*beta/((6.28001+beta)*nJacobi*nJacobi) ;
	root -= (fAbscissa[0]-root)*root1*root2*root3 ;
     }
     else if (i == nJacobi-1)
     {
	root1=(1.0+0.235002*beta)/(0.766001+0.118998*beta) ;
	root2=1.0/(1.0+0.639002*(nJacobi-4.0)/(1.0+0.71001*(nJacobi-4.0))) ;
	root3=1.0/(1.0+20.0*alpha/((7.5+alpha)*nJacobi*nJacobi)) ;
	root += (root-fAbscissa[nJacobi-4])*root1*root2*root3 ;
     } 
     else if (i == nJacobi) 
     {
	root1 = (1.0+0.37002*beta)/(1.67001+0.27998*beta) ;
	root2 = 1.0/(1.0+0.22*(nJacobi-8.0)/nJacobi) ;
	root3 = 1.0/(1.0+8.0*alpha/((6.28002+alpha)*nJacobi*nJacobi)) ;
	root += (root-fAbscissa[nJacobi-3])*root1*root2*root3 ;
     } 
     else
     {
	root = 3.0*fAbscissa[i-2]-3.0*fAbscissa[i-3]+fAbscissa[i-4] ;
     }
     alphaBeta = alpha + beta ;
     for (k=1;k<=maxNumber;k++)
     {
	temp = 2.0 + alphaBeta ;
	newton1 = (alpha-beta+temp*root)/2.0 ;
	newton2 = 1.0 ;
	for (j=2;j<=nJacobi;j++)
	{
	   newton3 = newton2 ;
	   newton2 = newton1 ;
	   temp = 2*j+alphaBeta ;
	   a = 2*j*(j+alphaBeta)*(temp-2.0) ;
 	   b = (temp-1.0)*(alpha*alpha-beta*beta+temp*(temp-2.0)*root) ;
	   c = 2.0*(j-1+alpha)*(j-1+beta)*temp ;
	   newton1 = (b*newton2-c*newton3)/a ;
	}
	newton = (nJacobi*(alpha - beta - temp*root)*newton1 +
	      2.0*(nJacobi + alpha)*(nJacobi + beta)*newton2)/
	     (temp*(1.0 - root*root)) ;
	rootTemp = root ;
	root = rootTemp - newton1/newton ;
	if (fabs(root-rootTemp) <= tolerance)
	{
	   break ;
	}
     }
     if (k > maxNumber) 
     {
        G4Exception("Too many iterations in G4GaussJacobiQ::G4GaussJacobiQ") ;
     }
     fAbscissa[i-1] = root ;
     fWeight[i-1] = exp(GammaLogarithm((G4double)(alpha+nJacobi)) + 
		        GammaLogarithm((G4double)(beta+nJacobi)) - 
		        GammaLogarithm((G4double)(nJacobi+1.0)) -
		        GammaLogarithm((G4double)(nJacobi + alphaBeta + 1.0)))
		        *temp*pow(2.0,alphaBeta)/(newton*newton2)                          ;
  }
}


// ----------------------------------------------------------
//
// Gauss-Jacobi method for integration of ((1-x)^alpha)*((1+x)^beta)*pFunction(x)
// from minus unit to plus unit .


G4double 
G4GaussJacobiQ::Integral() const 
{
   G4int i ;
   G4double integral = 0.0 ;
   for(i=0;i<fNumber;i++)
   {
      integral += fWeight[i]*fFunction(fAbscissa[i]) ;
   }
   return integral ;
}

