// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: testDataInterpolation.cc,v 1.3 1999-11-23 14:59:59 gcosmo Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
#include "G4ios.hh"
#include "globals.hh"
#include "G4DataInterpolation.hh"

G4double TestFunction(G4double x)
{
   return 10.0*exp(-0.1*x)*cos(x) ;
   //return 10.0*exp(-5*(x-pi)*(x-pi)) ;
}

int main()
{
   G4int i, j ;
   const G4int n = 15 ;
   G4double pY[n], pX[n], cof[n] ;
   G4double x, polcof, pol, yTest, deltaY, deltaX = twopi/n ;
   for(i=0;i<n;i++)
   {
      pX[i] = deltaX*i ;
      pY[i] = TestFunction(deltaX*i) ;
   }
   
   G4DataInterpolation myPolInt(pX,pY,n) ;
   myPolInt.PolIntCoefficient(cof) ;
   
   // Test PolCof against Polynom
   
   G4cout<<"Test function"<<"\t"<<"Delta Pol"<<"\t"<<"delta      "
       <<"\t"<<"Delta PolCof"<<"\t"<<"Delta Pol-PolCof"<<G4endl<<G4endl ;
   for(i=1;i<n-1;i++)
   {
      x = deltaX/2 + deltaX*i ;
      
      polcof = cof[n-1] ;
      for(j=n-2;j>=0;j--)
      {
	 polcof = polcof*x +cof[j] ;
      }
      pol = myPolInt.PolynomInterpolation(x,deltaY) ;
      G4cout<<TestFunction(x)<<"\t"
	  <<TestFunction(x) - pol<<"\t"
	  <<deltaY<<"\t"
	  <<TestFunction(x) - polcof<<"\t"
	  <<pol - polcof<<G4endl ;
   }
   G4cout<<G4endl ;
/* *************************
   
   // Test RationalPol against Polynomial
   
   G4cout<<"Test function"<<"\t"<<"Delta Pol"<<"\t"<<"delta      "
       <<"\t"<<"Delta RatPol"<<"\t"<<"delta"<<G4endl<<G4endl ;
   for(i=1;i<n-1;i++)
   {
      x = deltaX/2 + deltaX*i ;
      //yTest = myPolInt.RationalPolInterpolation(x,deltaY) ;
      G4cout<<TestFunction(x)<<"\t"
	  <<TestFunction(x) - myPolInt.PolynomInterpolation(x,deltaY)<<"\t"
	  <<deltaY<<"\t"
	  <<TestFunction(x) - myPolInt.RationalPolInterpolation(x,deltaY)<<"\t"
	  <<deltaY<<G4endl ;
   }
   // Test CubicSpline against Polynomial
   // Evaluation of start and finish first derivatives
   G4double deriStart = (pY[1]-pY[0])/deltaX ;
   G4double deriFinish = (pY[n-1]-pY[n-2])/deltaX ;
   
   G4DataInterpolation myPolIntCub(pX,pY,n,deriStart,deriFinish) ; // f''[i] is OK
   
   G4cout<<"Test function"<<"\t"<<"Delta Pol"<<"\t"<<"delta      "
       <<"\t"<<"Delta CubicSpline"<<"\t"<<"Delta FastCubicSpline"<<G4endl<<G4endl ;
   for(i=1;i<n-1;i++)
   {
      x = deltaX/2 + deltaX*i ;
      //yTest = myPolInt.RationalPolInterpolation(x,deltaY) ;
      G4cout<<TestFunction(x)<<"\t"
	  <<TestFunction(x) - myPolIntCub.PolynomInterpolation(x,deltaY)<<"\t"
	  <<deltaY<<"\t"
	  <<TestFunction(x) - myPolIntCub.CubicSplineInterpolation(x)<<"\t"
	  <<TestFunction(x) - myPolIntCub.FastCubicSpline(x,i)<<G4endl ;
   }
   G4cout<<G4endl ;
   G4cout<<"j"<<"\t"<<"x[j]"<<"\t"<<"pX"<<"Locate j"<<"\t"<<"Correlated j"<<G4endl ;
   G4int index ;
   for(i=1;i<n-1;i++)
   {
      x = deltaX/2 + deltaX*i ;
      index = i ;
      myPolInt.CorrelatedSearch(x,index) ;
      G4cout<<i<<"\t"<<pX[i]<<"\t"<<x<<"\t"
	  <<myPolInt.LocateArgument(x)<<"\t"<<index<<G4endl ;
   }
*/ ///////////////////////////   
   // myPolIntCub.~G4DataInterpolation() ;  
  return 0;
}
