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
// $Id: G4GaussLegendreQ.hh 67970 2013-03-13 10:10:06Z gcosmo $
//
// Class description:
//
// Class for Gauss-Legendre integration method
// Roots of ortogonal polynoms and corresponding weights are calculated based on
// iteration method (by bisection Newton algorithm). Constant values for initial
// approximations were derived from the book: M. Abramowitz, I. Stegun, Handbook
// of mathematical functions, DOVER Publications INC, New York 1965 ; chapters 9,
// 10, and 22 .
//
// ------------------------- CONSTRUCTORS: -------------------------------
//
// Constructor for GaussLegendre quadrature method. The value nLegendre set the
// accuracy required, i.e the number of points where the function pFunction will
// be evaluated during integration. The constructor creates the arrays for 
// abscissas and weights that used in Gauss-Legendre quadrature method. 
// The values a and b are the limits of integration of the pFunction.
// 
// G4GaussLegendreQ( function pFunction,
//		     G4int nLegendre           )
//
// -------------------------- METHODS:  ---------------------------------------
//
// Returns the integral of the function to be pointed by fFunction between a and b,
// by 2*fNumber point Gauss-Legendre integration: the function is evaluated exactly
// 2*fNumber Times at interior points in the range of integration. Since the weights
// and abscissas are, in this case, symmetric around the midpoint of the range of
// integration, there are actually only fNumber distinct values of each.
//
// G4double Integral(G4double a, G4double b) const 
//
// -----------------------------------------------------------------------
//
// Returns the integral of the function to be pointed by fFunction between a and b,
// by ten point Gauss-Legendre integration: the function is evaluated exactly
// ten Times at interior points in the range of integration. Since the weights
// and abscissas are, in this case, symmetric around the midpoint of the range of
// integration, there are actually only five distinct values of each
//
// G4double 
// QuickIntegral(G4double a, G4double b) const 
//
// ---------------------------------------------------------------------
//
// Returns the integral of the function to be pointed by fFunction between a and b,
// by 96 point Gauss-Legendre integration: the function is evaluated exactly
// ten Times at interior points in the range of integration. Since the weights
// and abscissas are, in this case, symmetric around the midpoint of the range of
// integration, there are actually only five distinct values of each
//
// G4double 
// AccurateIntegral(G4double a, G4double b) const 

// ------------------------------- HISTORY --------------------------------
//
// 13.05.97 V.Grichine (Vladimir.Grichine@cern.chz0

#ifndef G4GAUSSLEGENDREQ_HH
#define G4GAUSSLEGENDREQ_HH

#include "G4VGaussianQuadrature.hh"

class G4GaussLegendreQ : public G4VGaussianQuadrature
{
public:
        explicit G4GaussLegendreQ( function pFunction ) ;
        

        G4GaussLegendreQ( function pFunction,
			  G4int nLegendre           ) ;
			       
        // Methods
			     
        G4double Integral(G4double a, G4double b) const ;

        G4double QuickIntegral(G4double a, G4double b) const ;
                              
        G4double AccurateIntegral(G4double a, G4double b) const ;

private:

	G4GaussLegendreQ(const G4GaussLegendreQ&);
	G4GaussLegendreQ& operator=(const G4GaussLegendreQ&);
};

#endif
