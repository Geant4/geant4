//  Bogacki-Shampine - 4 - 3(2) non-FSAL implementation 
//
//  An implementation of the embedded RK method from the paper 
// [1] P. Bogacki and L. F. Shampine, “A 3(2) pair of Runge - Kutta formulas,” 
// Appl. Math. Lett., vol. 2, no. 4, pp. 321–325, Jan. 1989.
//
//  This version does not utilise the FSAL property of the method,
//  which would allow the reuse of the last derivative in the next step.
//  (Alternative FSAL implementation created with revised interface)
//
//  Implemented by Somnath Banerjee
// Work supported by the Google Summer of Code 2015.
//  Supervision / code review: John Apostolakis
//
// First version: 20 May 2015
//
// This code is made available subject to the Geant4 license, a copy of
// which is available at
//   http://geant4.org/license
//
///////////////////////////////////////////////////////////////////////////////


#ifndef G4BOGACKI_SHAMPINE23_H
#define G4BOGACKI_SHAMPINE23_H

#include "G4MagIntegratorStepper.hh"

class G4BogackiShampine23 : public G4MagIntegratorStepper{


 public:
 	//constructor
 	G4BogackiShampine23( G4EquationOfMotion *EqRhs,
                     G4int numberOfVariables = 6,
                     G4bool primary= true ) ;

 	//destructor
 	~G4BogackiShampine23() ;

 	//Stepper
 	 void Stepper(const G4double y[],
                  const G4double dydx[],
                        G4double h,
                        G4double yout[],
                        G4double yerr[] ) ;

    G4double  DistChord()   const;
    G4int IntegratorOrder() const { return 2; }
     G4bool isFSAL() const{ return true; }
	G4double *getLastDydx();
    
	G4BogackiShampine23(const G4BogackiShampine23&);
   G4BogackiShampine23& operator=(const G4BogackiShampine23&);


 private:

    G4double *ak2, *ak3, *ak4, *yTemp, *yIn;
      // for storing intermediate 'k' values in stepper
    G4double *pseudoDydx_for_DistChord;
    
    G4double fLastStepLength;
    G4double *fLastInitialVector, *fLastFinalVector,
             *fLastDyDx, *fMidVector, *fMidError;
      // for DistChord calculations

    G4BogackiShampine23* fAuxStepper;

//	G4int No_of_vars;
//	G4double hinit, tinit, tmax, *yinit;
//	double hmax, hmin, safe_const, err0, Step_factor;
//	void (*derivs)(double, double *, double *);

};

#endif /* G4BogackiShampine23 */





























