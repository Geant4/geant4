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
// -------------------------------------------------------------------
//      GEANT4 Class file
//
//
//      File name:     G4Clebsch
//
//      Author:        Maria Grazia Pia (MariaGrazia.Pia@genova.infn.it)
// 
//      Creation date: 15 April 1999
//
//      Modifications: February 2015, Jason Detwiler (jasondet@gmail.com)
//                     - Update variable names to clarify that twice 
//                       the spin value is being used.
//                     - Add function to obtain raw CG coeffients
//                     - Speed up Wigner 3j evaluation by removing 
//                       reallocation of vector<double>'s on every call
//                     - Use G4Pow instead of internal logs variable
//                     - Add functions for Wigner 6j and 9j symbols, 
//                       Racah coefficients, and d-Matrices
//                     - Make functions statically available
//                     - Eliminate unnecessary constructors, destructors,
//                       operators (class has no members, so default versions
//                       are adequate)
//                     11.06.2015 V.Ivanchenko adopted for Geant4 source 
//      
//      Clebsch-Gordan coefficient and related algebra
//
// -------------------------------------------------------------------

#ifndef G4CLEBSCH_HH
#define G4CLEBSCH_HH

#include "globals.hh"
#include <vector>

class G4Clebsch 
{
public:
  // Calculates the standard Clebsch-Gordan coefficient with Condon-Shortley
  // phase convention.
  static 
  G4double ClebschGordanCoeff(G4int twoJ1, G4int twoM1, 
			      G4int twoJ2, G4int twoM2, 
			      G4int twoJ);

  // Calculates the square of the corresponding Clebsch-Gordan coefficient 
  static 
  G4double ClebschGordan(G4int twoJ1, G4int twoM1, 
			 G4int twoJ2, G4int twoM2, 
			 G4int twoJ);

  static 
  std::vector<G4double> GenerateIso3(G4int twoJ1, G4int twoM1, 
				     G4int twoJ2, G4int twoM2, 
				     G4int twoJOut1, G4int twoJOut2);

  static 
  G4double Weight(G4int twoJ1, G4int twoM1, 
		  G4int twoJ2, G4int twoM2, 
		  G4int twoJOut1, G4int twoJOut2);

  // Deprecated G4double version of Wigner3J that assumes user sent in half-int
  // values. Use G4int version instead, which guarantees validity of arguments.
  static 
  G4double Wigner3J(G4double j1, G4double j2, G4double j3, 
		    G4double m1, G4double m2, G4double m3);
  
  // Wigner's 3J symbols (CG-coeffs with a different phase and normalization)
  static G4double Wigner3J(G4int twoJ1, G4int twoM1, 
                           G4int twoJ2, G4int twoM2, 
                           G4int twoJ3);
  static G4double Wigner3J(G4int twoJ1, G4int twoM1, 
                           G4int twoJ2, G4int twoM2, 
                           G4int twoJ3, G4int twoM3);

  // Calculates the normalized Clebsch-Gordan coefficient, that is the prob 
  // of isospin decomposition of (J,m) into J1, J2, m1, m2
  static 
  G4double NormalizedClebschGordan(G4int twoJ,  G4int twom, 
				   G4int twoJ1, G4int twoJ2,
				   G4int twom1, G4int twom2);
  
  // Triangle coefficient appearing in explicit calculation of CG-coeffs
  static G4double TriangleCoeff(G4int twoA, G4int twoB, G4int twoC);

  // Wigner's 6J symbols (generalization of CG's for coupling of 3 angular momenta). 
  static G4double Wigner6J(G4int twoJ1, G4int twoJ2, G4int twoJ3, 
		           G4int twoJ4, G4int twoJ5, G4int twoJ6);

  // Racah's W Coeffs (6J's with a different phase and normalization)
  static G4double RacahWCoeff(G4int twoJ1,  G4int twoJ2, 
		              G4int twoJ,   G4int twoJ3, 
		              G4int twoJ12, G4int twoJ23);

  // Wigner's 9J symbols (generalization of CG's for coupling of 4 angular momenta). 
  static G4double Wigner9J(G4int twoJ1, G4int twoJ2, G4int twoJ3, 
		           G4int twoJ4, G4int twoJ5, G4int twoJ6,
		           G4int twoJ7, G4int twoJ8, G4int twoJ9);

  // Wigner's little-d matrix
  static G4double WignerLittleD(G4int twoJ, G4int twoM, G4int twoN, 
				G4double cosTheta);
};

#endif

