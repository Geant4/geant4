// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: testIntegration.cc,v 1.2 1999-11-16 17:31:21 gcosmo Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// Test program for G4SimpleIntegration class. The function exp(-x)*cos(x) is
// integrated between zero and two pi. The true result is 0.499066278634
//

#include "G4ios.hh"
#include "globals.hh"
#include "G4SimpleIntegration.hh"

G4double TestFunction(G4double x)
{
  return exp(-x)*cos(x) ;
}

int main()
{
   G4int i, n ;
   G4double pTolerance ;
   G4double a = 0.0 ;
   G4double b = twopi ;
   G4SimpleIntegration myIntegrand(TestFunction) ;

   G4cout<<"Iteration"<<"\t"<<"Trapezoidal"<<"\t"
       <<"MidPoint"<<"\t"<<"Gauss"<<"\t"<<"Simpson"<<endl ;
   for(i=0;i<13;i++)
   {
      n = (int)pow(2,i) ;
      G4cout<<n<<"\t"
	  <<myIntegrand.Trapezoidal(a,b,n)<<"\t"
	  <<myIntegrand.MidPoint(a,b,n)<<"\t"
	  <<myIntegrand.Gauss(a,b,n)<<"\t"
	  <<myIntegrand.Simpson(a,b,n)<<endl ;
   }
   G4cout<<endl ;
   for(i=0;i<13;i++)
   {
      pTolerance = pow(10.0,-i) ;
      G4SimpleIntegration adaptIntegrand(TestFunction,pTolerance) ;
      G4cout<<adaptIntegrand.AdaptGaussIntegration(a,b)<<endl;
   }
   return 0;
}
