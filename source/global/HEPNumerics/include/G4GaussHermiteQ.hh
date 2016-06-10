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
// $Id: G4GaussHermiteQ.hh 67970 2013-03-13 10:10:06Z gcosmo $
//
// Class description:
//
// Roots of ortogonal polynoms and corresponding weights are calculated based on
// iteration method (by bisection Newton algorithm). Constant values for initial
// approximations were derived from the book: M. Abramowitz, I. Stegun, Handbook
// of mathematical functions, DOVER Publications INC, New York 1965 ; chapters 9,
// 10, and 22 .
//
// --------------------------------------------------------------------------
//
// Constructor for Gauss-Hermite quadrature method . The function GaussHermite
// should be called then
//
// G4GaussHermiteQ( function pFunction, G4int nHermite  ) 
//
// ----------------------------------------------------------------------------
//
// Gauss-Hermite method for integration of std::exp(-x*x)*nFunction(x) from minus infinity
// to plus infinity .
//
// G4double Integral() const 

// ------------------------------- HISTORY -------------------------------------
//
// 13.05.97 V.Grichine (Vladimir.Grichine@cern.chz0

#ifndef G4GAUSSHERMITEQ_HH
#define G4GAUSSHERMITEQ_HH

#include "G4VGaussianQuadrature.hh"

class G4GaussHermiteQ : public G4VGaussianQuadrature
{
public:
        // Constructor

        G4GaussHermiteQ( function pFunction, G4int nHermite  ) ;
			       
        // Methods
			     
        G4double Integral() const ;


private:

	G4GaussHermiteQ(const G4GaussHermiteQ&);
	G4GaussHermiteQ& operator=(const G4GaussHermiteQ&);

};

#endif
