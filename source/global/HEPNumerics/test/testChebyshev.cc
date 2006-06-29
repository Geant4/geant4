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
// $Id: testChebyshev.cc,v 1.6 2006-06-29 19:00:28 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// Test program for G4ChebyshevApproximation class. The function std::exp(-x)*std::cos(x) is
// integrated between zero and two pi. The true result is 0.499066278634
//

#include "G4ios.hh"
#include "globals.hh"
#include "G4ChebyshevApproximation.hh"

G4double TestChebyshev(G4double x)
{
  return std::sqrt(1-x*x)*std::cos(x) ;
}

G4double TestFunction(G4double x)
{
  return std::exp(-x)*std::cos(x) ;
}

G4double TestHermite(G4double x)
{
  return x*x*std::cos(x) ;
}

G4double ExpFunction(G4double x)
{
  return std::exp(x) ;
}

G4double SinFunction(G4double x)
{
  return std::sin(x) ;
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
	  <<test<<"\t"<<tolerance<<G4endl ;
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
	  <<test<<"\t"<<tolerance<<G4endl ;
   }
   G4double test, tolerance, true =  TestFunction(x) ;
   for(i=1;i<n;i++)
   {
      G4ChebyshevApproximation myChebyshev(TestFunction,i,a,b) ;
      test = myChebyshev.ChebyshevEvaluation(x) ;
      tolerance = 2*(true-test)/(true+test) ;
      G4cout<<"n = "<<i<<"\t"<<"true = "<<true
	  <<" and n-point ChebEval =  "
	  <<test<<"\t"<<tolerance<<G4endl ;
   }
*/ ///////////////////////////////////////////////////   
  return 0;
}
