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
// $Id: G4AnalyticalPolSolver.hh,v 1.4 2005/05/19 07:37:10 gcosmo Exp $
// GEANT4 tag $Name: geant4-08-00 $
//
// Class description:
//
// G4AnalyticalPolSolver allows the user to solve analytically a polynomial
// equation up to the 4th order. This is used by CSG solid tracking functions
// like G4Torus.
//
// The algorithm has been adapted from the CACM Algorithm 326:
//
//   Roots of low order polynomials
//   Author: Terence R.F.Nonweiler
//   CACM  (Apr 1968) p269
//   Translated into C and programmed by M.Dow
//   ANUSF, Australian National University, Canberra, Australia
//   m.dow@anu.edu.au
//
// Suite of procedures for finding the (complex) roots of the quadratic,
// cubic or quartic polynomials by explicit algebraic methods.
// Each Returns:
//
//   x=r[1][k] + i r[2][k]  k=1,...,n, where n={2,3,4}
//
// as roots of:
// sum_{k=0:n} p[k] x^(n-k) = 0
// Assumes p[0] != 0. (< or > 0) (overflows otherwise)

// --------------------------- HISTORY --------------------------------------
//
// 13.05.05 V.Grichine ( Vladimir.Grichine@cern.ch )
//          First implementation in C++

#ifndef G4AN_POL_SOLVER_HH
#define G4AN_POL_SOLVER_HH

#include  "G4Types.hh"

class G4AnalyticalPolSolver 
{
  public:  // with description

    G4AnalyticalPolSolver();
    ~G4AnalyticalPolSolver();

    G4int QuadRoots(    G4double p[5], G4double r[3][5]);
    G4int CubicRoots(   G4double p[5], G4double r[3][5]);
    G4int BiquadRoots(  G4double p[5], G4double r[3][5]);
    G4int QuarticRoots( G4double p[5], G4double r[3][5]);
};

#endif
