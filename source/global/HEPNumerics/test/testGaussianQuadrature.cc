// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: testGaussianQuadrature.cc,v 1.3 1999-11-23 15:00:00 gcosmo Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// Test program for G4GaussianQuadrature class. The function exp(-x)*cos(x) is
// integrated between zero and two pi. The true result is 0.499066278634
//

#include "G4ios.hh"
#include "globals.hh"
// #include "G4VGaussianQuadrature.hh"
#include "G4GaussChebyshevQ.hh"
#include "G4GaussHermiteQ.hh"
#include "G4GaussJacobiQ.hh"
#include "G4GaussLaguerreQ.hh"
#include "G4GaussLegendreQ.hh"

G4double TestChebyshev(G4double x)
{
  return sqrt(1-x*x)*cos(x) ;
}

G4double TestFunction(G4double x)
{
  return exp(-x)*cos(x) ;
}

G4double TestHermite(G4double x)
{
  return x*x*cos(x) ;
}

G4double CosFunction(G4double x)
{
  return cos(x) ;
}

int main()
{
    G4int i, n = 20;
   // G4double pTolerance ;
   G4double a = -1.0 ;
   G4double b = 1.0 ;
   // G4double true =  pi*0.4400505857 ;
   for(i=1;i<20;i++)
   {
      n = 1*i ;
      G4GaussChebyshevQ myChebyshev(TestChebyshev,n) ;
      G4cout<<"n = "<<n<<"\t"<<"true = "<<true<<"  and n-point Gauss-Chebyshev =  "
	  <<myChebyshev.Integral(a,b)<<G4endl ;
   }
   // G4double true = pi*0.4400505857 ;
   for(i=1;i<20;i++)
   {
      n = 1*i ;
      G4GaussJacobiQ myJacobi(CosFunction,0.5,0.5,n) ;
      G4cout<<"n = "<<n<<"\t"<<"true = "<<true<<"  and n-point Gauss-Jacobi =  "
	  <<myJacobi.Integral()<<G4endl ;
   }
   G4double true = 2*0.125*sqrt(pi)*exp(-0.25) ;
   for(i=1;i<20;i++)
   {
      n = 1*i ;
      G4GaussHermiteQ myHermite(TestHermite,n) ;
      G4cout<<"n = "<<n<<"\t"<<"true = "<<true<<"  and n-point GaussHermite =  "
	  <<myHermite.Integral()<<G4endl ;
   }
   /* *******************
   G4GaussianQuadrature myHermite(n) ;
   for(i=0;i<(n+1)/2;i++)
   {
      G4cout<<i<<"\t"<<myHermite.GetAbscissa(i)<<"\t"
	  <<myHermite.GetWeight(i)<<G4endl ;
   }
   G4GaussianQuadrature myLaguerre(0.0,n) ;
   for(i=0;i<n;i++)
   {
      G4cout<<i<<"\t"<<myLaguerre.GetAbscissa(i)<<"\t"
	  <<myLaguerre.GetWeight(i)<<G4endl ;
   }
   */ /////////////////////////////////////
   for(i=1;i<20;i++)
   {
      n = 1*i ;
      G4GaussLaguerreQ myLaguerre(CosFunction,0.0,n) ;
      G4cout<<"n = "<<n<<"\t"<<"true = 0.5 "<<"  and n-point GaussLaguerre =  "
	  <<myLaguerre.Integral()<<G4endl ;
   }
   G4GaussLegendreQ myIntegrand(TestFunction) ;
   G4cout<<"true is 0.499066278634 "<<"  and QuickGaussLegendre is  "<<
      myIntegrand.QuickIntegral(0,2*pi)<<G4endl ;
   G4cout<<"true is 0.499066278634 "<<"  and AccurateGaussLegendre is  "<<
      myIntegrand.AccurateIntegral(0,2*pi)<<G4endl ;


   for(i=1;i<20;i++)
   {
      n = 8*i ;
      G4GaussLegendreQ myLegendre(TestFunction, n) ;
      G4cout<<myLegendre.GetNumber()<<
      "true is 0.5 "<<"  and n-point GaussLegendre is  "
	  <<myLegendre.Integral(0,200*pi)<<G4endl ;
   }
   /* **************************************************
   
   */ /////////////////////////////////////
  return 0;
}
