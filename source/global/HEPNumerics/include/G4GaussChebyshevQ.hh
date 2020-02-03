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
//
// Class description:
//
// Class for Gauss-Chebyshev quadrature method
// Roots of ortogonal polynoms and corresponding weights are calculated based on
// iteration method (by bisection Newton algorithm). Constant values for initial
// approximations were derived from the book: M. Abramowitz, I. Stegun, Handbook
// of mathematical functions, DOVER Publications INC, New York 1965 ; chapters 9,
// 10, and 22 .
//
// ------------------------------ CONSTRUCTORS ----------------------------
//
// Constructor for Gauss-Chebyshev quadrature method
//
// G4GaussChebyshevQuadrature( function pFunction,
//			       G4int nChebyshev       )
//
//
//
// ------------------------------- METHODS -----------------------------------
//
// Integrates function pointed by fFunction from a to b by Gauss-Chebyshev quadrature
// method
//
// G4double Integral(G4double a, G4double b) const 

// ------------------------------- HISTORY --------------------------------
//
// 13.05.97 V.Grichine (Vladimir.Grichine@cern.ch)

#ifndef G4GAUSSCHEBYSHEVQ_HH
#define G4GAUSSCHEBYSHEVQ_HH

#include "G4VGaussianQuadrature.hh"

class G4GaussChebyshevQ : public G4VGaussianQuadrature
{
public:
        // Constructor/destructor

        G4GaussChebyshevQ( function pFunction,
			   G4int nChebyshev           ) ;
			   
	~G4GaussChebyshevQ() ;		   
			       
        // Methods
			     
        G4double Integral(G4double a, G4double b) const ;


private:

	G4GaussChebyshevQ(const G4GaussChebyshevQ&);
	G4GaussChebyshevQ& operator=(const G4GaussChebyshevQ&);

};

#endif
