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
// $Id: G4GaussJacobiQ.hh,v 1.4 2001-07-11 10:00:39 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
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
