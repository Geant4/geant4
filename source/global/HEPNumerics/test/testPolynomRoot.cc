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
// Test program for G4JTPolynomialSolver class. 
//

#include "G4ios.hh"
#include "globals.hh"
#include "G4JTPolynomialSolver.hh"

int main()
{

  G4cout << "root finding starts ... " << G4endl ;

  //  double c[5] = { 1, 0.99778585, 26.793624, 17.923849, -1.6282041} ;
  //  double c[5] = { 1, 1.3062095 , -1.8130575, -8.5895246 , -2.5500752  } ;
  //  double c[5] = {1, -0.56733739, 1.8378874 , 3.9504828 , 0.61193337 } ;
  //  double c[5] = { 1,-3.1634943 , -1.1041289, 10.876123 ,-3.6351417 } ;
  //  double c[5] = { 1, 7.5099001 , 26.263224, -14.43852 , -13.015633 } ;
  //  double c[5] = { 1, -0.53770864, -7.7970002, -4.4380489 ,-0.65266367 } ;

  G4double c[5] = { 1, 5.0977493 , 3.0361937, -11.452036, -2.6873796 } ;
  G4double sr[4] ;  // real part
  G4double si[4] ;  // imaginary part

  G4int degree = 4 ;
  G4int num ;       // number of roots 

  // solve the polynom analytically
  G4JTPolynomialSolver trapEq ;
  num = trapEq.FindRoots(c, degree, sr, si) ;
  

  G4cout.precision(16) ;
  G4cout << "number of roots found = " << num << G4endl ;

  for ( int i = 0 ; i < num ; i++ ) {
    G4cout << "solution " << i << " = " << sr[i] << " + " << si[i] << "i" << G4endl ;
  }
  

  return 1 ;

}
