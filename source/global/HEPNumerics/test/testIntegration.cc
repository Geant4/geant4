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
// $Id: testIntegration.cc,v 1.6 2006-06-29 19:00:39 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// Test program for G4SimpleIntegration class. The function std::exp(-x)*std::cos(x) is
// integrated between zero and two pi. The true result is 0.499066278634
//

#include "G4ios.hh"
#include "globals.hh"
#include "G4SimpleIntegration.hh"

G4double TestFunction(G4double x)
{
  return std::exp(-x)*std::cos(x) ;
}

int main()
{
   G4int i, n ;
   G4double pTolerance ;
   G4double a = 0.0 ;
   G4double b = twopi ;
   G4SimpleIntegration myIntegrand(TestFunction) ;

   G4cout<<"Iteration"<<"\t"<<"Trapezoidal"<<"\t"
       <<"MidPoint"<<"\t"<<"Gauss"<<"\t"<<"Simpson"<<G4endl ;
   for(i=0;i<13;i++)
   {
      n = (int)pow(2,i) ;
      G4cout<<n<<"\t"
	  <<myIntegrand.Trapezoidal(a,b,n)<<"\t"
	  <<myIntegrand.MidPoint(a,b,n)<<"\t"
	  <<myIntegrand.Gauss(a,b,n)<<"\t"
	  <<myIntegrand.Simpson(a,b,n)<<G4endl ;
   }
   G4cout<<G4endl ;
   for(i=0;i<13;i++)
   {
      pTolerance = std::pow(10.0,-i) ;
      G4SimpleIntegration adaptIntegrand(TestFunction,pTolerance) ;
      G4cout<<adaptIntegrand.AdaptGaussIntegration(a,b)<<G4endl;
   }
   return 0;
}
