//
// ********************************************************************
// * DISCLAIMER                                                       *
// *                                                                  *
// * The following disclaimer summarizes all the specific disclaimers *
// * of contributors to this software. The specific disclaimers,which *
// * govern, are listed with their locations in:                      *
// *   http://cern.ch/geant4/license                                  *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.                                                             *
// *                                                                  *
// * This  code  implementation is the  intellectual property  of the *
// * GEANT4 collaboration.                                            *
// * By copying,  distributing  or modifying the Program (or any work *
// * based  on  the Program)  you indicate  your  acceptance of  this *
// * statement, and all its terms.                                    *
// ********************************************************************
//
//
// Test program for G4PolynomRoot class. 
//

#include "G4ios.hh"
#include "globals.hh"
#include "G4PolynomRoot.hh"

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
  G4PolynomRoot trapEq ;
  num = trapEq.FindRoots(c, degree, sr, si) ;
  

  G4cout.precision(16) ;
  G4cout << "number of roots found = " << num << G4endl ;

  for ( int i = 0 ; i < num ; i++ ) {
    G4cout << "solution " << i << " = " << sr[i] << " + " << si[i] << "i" << G4endl ;
  }
  

  return 1 ;

}
