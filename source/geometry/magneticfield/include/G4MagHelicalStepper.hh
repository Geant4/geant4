// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4MagHelicalStepper.hh,v 1.3 2000-04-12 18:28:51 japost Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//                                         Started from G4MagErrorStepper.hh
//
// Abstract base class (ie Interface)
// -------------------
//    for integrator of particle's equation of motion,
//    used in tracking in space dependent magnetic field
// -----------------------------------------------------
//
// History:
// 05.11.98  J.Apostolakis   Creation of new ABC 
//
#ifndef G4MagHelicalStepper_hh
#define G4MagHelicalStepper_hh
#include "globals.hh"
#include "G4MagIntegratorStepper.hh"
#include "G4Mag_EqRhs.hh"
#include "G4ThreeVector.hh"

class G4MagHelicalStepper : public G4MagIntegratorStepper
{
 public:

  G4MagHelicalStepper(G4Mag_EqRhs *EqRhs);
  ~G4MagHelicalStepper(){} ;

  //   The stepper for the Runge Kutta integration. The stepsize 
  //     is fixed, equal to h.
  //     Integrates ODE starting values y[0 to 6 ]
  //     Outputs yout[] and its estimated error yerr[].
  
  void  Stepper(  const G4double y[],
		  const G4double dydx[],
		  const G4double h,
		  G4double yout[],
		  G4double yerr[]  );

  // performs a 'dump' Step without error calculation.
  
  virtual  void  DumbStepper(  const G4double y[],
			       G4ThreeVector   Bfld,
			       G4double  h,
			       G4double yout[] ) = 0;

  // Estimate maximum distance of curved solution and chord ... 
  
  G4double DistChord()   const;

 //  --- Methods used to implement all the derived classes -----
 protected:

  // a linear Step in regions without magnetic field

  inline void LinearStep( const G4double  yIn[],
			  const G4double  h,
			  G4double  yHelix[]);
  
  // a first order Step along a helix inside the field

  void AdvanceHelix( const G4double  yIn[],
		     G4ThreeVector   Bfld,
		     G4double  h,
		     G4double  yHelix[]);    // output 

  // evaluate the field at a certain point

  // void MagFieldEvaluate( const G4double y[], G4double B[] )   
  //  { GetEquationOfMotion()->  GetFieldValue(y, B); }

  void MagFieldEvaluate( const G4double y[], G4ThreeVector& Bfield )   
    { G4double B[3];
      GetEquationOfMotion()->  GetFieldValue(y, B);
      Bfield= G4ThreeVector( B[0], B[1], B[2] );
    }

 private:
  
  // Data stored in order to find the chord
  G4ThreeVector yInitial, yMidPoint, yFinal;

  G4Mag_EqRhs*  fPtrMagEqOfMot;
};

#include  "G4MagHelicalStepper.icc"

#endif  /* G4MagHelicalStepper_hh */
