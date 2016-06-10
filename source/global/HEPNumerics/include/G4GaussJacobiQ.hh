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
// $Id: G4GaussJacobiQ.hh 67970 2013-03-13 10:10:06Z gcosmo $
//
// Class description:
//
// Roots of ortogonal polynoms and corresponding weights are calculated based on
// iteration method (by bisection Newton algorithm). Constant values for initial
// approximations were derived from the book: M. Abramowitz, I. Stegun, Handbook
// of mathematical functions, DOVER Publications INC, New York 1965 ; chapters 9,
// 10, and 22 .
//
// ---------------------------------------------------------------------------
//
// Constructor for Gauss-Jacobi integration method. 
//
// G4GaussJacobiQ( function pFunction,
//                 G4double alpha,
//                 G4double beta, 
//		   G4int nJacobi   ) 
//
// ----------------------------------------------------------------------------
//
// Gauss-Jacobi method for integration of ((1-x)^alpha)*((1+x)^beta)*pFunction(x)
// from minus unit to plus unit .
//
// G4double Integral() const

// ------------------------------- HISTORY -------------------------------------
//
// 13.05.97 V.Grichine (Vladimir.Grichine@cern.chz0

#ifndef G4GAUSSJACOBIQ_HH
#define G4GAUSSJACOBIQ_HH

#include "G4VGaussianQuadrature.hh"

class G4GaussJacobiQ : public G4VGaussianQuadrature
{
public:
        // Constructor

        G4GaussJacobiQ( function pFunction, 
	                G4double alpha,
	                G4double beta,
	                G4int nJacobi         ) ;
			       
        // Methods
			     
        G4double Integral() const ;

private:

	G4GaussJacobiQ(const G4GaussJacobiQ&);
	G4GaussJacobiQ& operator=(const G4GaussJacobiQ&);
};

#endif
