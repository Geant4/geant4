// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: testChebyshev.cc,v 1.1 1999-01-07 16:08:57 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// Test program for G4ChebyshevApproximation class. The function exp(-x)*cos(x) is
// integrated between zero and two pi. The true result is 0.499066278634
//

#include "G4ios.hh"
#include "globals.hh"
#include "G4ChebyshevApproximation.hh"

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

G4double ExpFunction(G4double x)
{
  return exp(x) ;
}

G4double SinFunction(G4double x)
{
  return sin(x) ;
}

main()
{
   G4int i, k, m, n = 30;
   G4double x = 3.0 ;
   G4double a = 0.0 ;
   G4double b = 10.0 ;
   
   G4double test, tolerance, true =  ExpFunction(x) - 1.0 ;
   for(i=5;i<=n;i++)
   {
      G4ChebyshevApproximation myChebyshev(ExpFunction,a,b,i) ;  // integral
      test = myChebyshev.ChebyshevEvaluation(x) ;
      tolerance = 2*(true-test)/(true+test) ;
      G4cout<<"n = "<<i<<"\t"<<"true = "<<true
	  <<" and n-point ChebEval =  "
	  <<test<<"\t"<<tolerance<<endl ;
   }
   /* ********************************************************
   n = k*m + 10 ;
   
   G4double x = 3.0 ;
   G4double delta = m*0.1*x ;
   G4double a = x - delta ;
   G4double b = x + delta ;
   
   G4double test, tolerance, true =  SinFunction(x) ;
   for(i=1+m;i<=n;i++)
   {
      G4ChebyshevApproximation myChebyshev(SinFunction,i,m,a,b) ; // m-derivative
      test = myChebyshev.ChebyshevEvaluation(x) ;
      tolerance = 2*(true-test)/(true+test) ;
      G4cout<<"n = "<<i<<"\t"<<"true = "<<true
	  <<" and n-point ChebEval =  "
	  <<test<<"\t"<<tolerance<<endl ;
   }
   G4double test, tolerance, true =  TestFunction(x) ;
   for(i=1;i<n;i++)
   {
      G4ChebyshevApproximation myChebyshev(TestFunction,i,a,b) ;
      test = myChebyshev.ChebyshevEvaluation(x) ;
      tolerance = 2*(true-test)/(true+test) ;
      G4cout<<"n = "<<i<<"\t"<<"true = "<<true
	  <<" and n-point ChebEval =  "
	  <<test<<"\t"<<tolerance<<endl ;
   }
*/ ///////////////////////////////////////////////////   
  return 0;
}
