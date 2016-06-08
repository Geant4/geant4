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
// $Id: G4Solver.hh,v 1.3.2.1 2001/06/28 19:13:04 gunter Exp $
// GEANT4 tag $Name:  $
//
// Hadronic Process: Nuclear De-excitations
// by V. Lara

#ifndef G4Solver_h
#define G4Solver_h 1

#include "globals.hh"

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
// 	template <class Function> G4bool Bisection(const Function & theFunction);
    G4bool Bisection(Function & theFunction);	
	
    // Calculates the root by the Regula-Falsi method
// 	template <class Function> G4bool RegulaFalsi(const Function & theFunction);
    G4bool RegulaFalsi(Function & theFunction);
	
	
    // Calculates the root by the Brent's method
    G4bool Brent(Function & theFunction);
	
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
