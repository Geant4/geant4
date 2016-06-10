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
// $Id: G4AnalyticalPolSolver.hh 67970 2013-03-13 10:10:06Z gcosmo $
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
