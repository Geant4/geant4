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
// $Id: G4GaussLegendreQ.hh,v 1.4 2001-07-11 10:00:39 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
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
        G4GaussLegendreQ( function pFunction ) ;
        

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
