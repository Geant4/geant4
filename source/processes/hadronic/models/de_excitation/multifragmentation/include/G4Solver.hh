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
// $Id: G4Solver.hh 67983 2013-03-13 10:42:03Z gcosmo $
//
// Hadronic Process: Nuclear De-excitations
// by V. Lara

#ifndef G4Solver_h
#define G4Solver_h 1

#include "globals.hh"

#include <cmath>

#define DefaultTolerance 5.0e-14

template <class Function> class G4Solver 
{
public:
    enum {DefaultMaxIter = 100};
	
    // default constructor
    G4Solver() : MaxIter(DefaultMaxIter), tolerance(DefaultTolerance),
		 a(0.0), b(0.0), root(0.0) {};
	
    G4Solver(const G4int iterations, const G4double tol) :
	MaxIter(iterations), tolerance(tol),
	a(0.0), b(0.0), root(0.0) {};

    // copy constructor	
    G4Solver(const G4Solver & right);

    // destructor
    ~G4Solver() {};
	
    // operators
    G4Solver & operator=(const G4Solver & right);
    G4bool operator==(const G4Solver & right) const;
    G4bool operator!=(const G4Solver & right) const;
		
    G4int GetMaxIterations(void) const {return MaxIter;}
    void SetMaxIterations(const G4int iterations) {MaxIter=iterations;}
	
    G4double GetTolerance(void) const {return tolerance;}
    void SetTolerance(const G4double epsilon) {tolerance = epsilon;}
	
	
    G4double GetIntervalLowerLimit(void) const {return a;}
    G4double GetIntervalUpperLimit(void) const {return b;}
	
    void SetIntervalLimits(const G4double Limit1, const G4double Limit2);
    
    G4double GetRoot(void) const {return root;}
    
    // Calculates the root by the Bisection method
    G4bool Bisection(Function & theFunction);	
	
    // Calculates the root by the Regula-Falsi method
    G4bool RegulaFalsi(Function & theFunction);
	
	
    // Calculates the root by the Brent's method
    G4bool Brent(Function & theFunction);

    // Calculates the root by the Inverse Parabolic Interpolation method 
    // due to Jack Crenshaw
    G4bool Crenshaw(Function & theFunction);
	
private:

    // Maximum number of iterations
    G4int MaxIter;

    // 
    G4double tolerance;

    // interval limits [a,b] which should bracket the root
    G4double a;
    G4double b;

    // The root
    G4double root;

};

#include "G4Solver.icc"

#endif
