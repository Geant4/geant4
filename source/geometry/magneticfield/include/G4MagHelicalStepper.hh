// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4MagHelicalStepper.hh,v 1.6 2001-03-22 18:48:52 japost Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
//
// class G4MagHelicalStepper
//
// Class description:
//
// Abstract base class for integrator of particle's equation of motion,
// used in tracking in space dependent magnetic field

// History:
// - 05.11.98  J.Apostolakis   Creation of new ABC 

#ifndef G4MagHelicalStepper_hh
#define G4MagHelicalStepper_hh

#include "globals.hh"
#include "G4MagIntegratorStepper.hh"
#include "G4Mag_EqRhs.hh"
#include "G4ThreeVector.hh"

class G4MagHelicalStepper : public G4MagIntegratorStepper
{
  public:  // with description

    G4MagHelicalStepper(G4Mag_EqRhs *EqRhs);
    virtual ~G4MagHelicalStepper();
  
    void Stepper( const G4double y[],
		  const G4double dydx[],
		        G4double h,
		        G4double yout[],
		        G4double yerr[]  );
      // The stepper for the Runge Kutta integration.
      // The stepsize is fixed, equal to h.
      // Integrates ODE starting values y[0 to 6]
      // Outputs yout[] and its estimated error yerr[].
  
    virtual  void DumbStepper( const G4double y[],
			       G4ThreeVector   Bfld,
			       G4double  h,
			       G4double yout[] ) = 0;
      // Performs a 'dump' Step without error calculation.
  
    G4double DistChord() const;
      // Estimate maximum distance of curved solution and chord ... 

  protected:  // with description

    inline void LinearStep( const G4double  yIn[],
	                          G4double  h,
			          G4double  yHelix[]);
      // A linear Step in regions without magnetic field.

    void AdvanceHelix( const G4double  yIn[],
		       G4ThreeVector   Bfld,
		       G4double  h,
		       G4double  yHelix[]);    // output 
      // A first order Step along a helix inside the field.

    inline void MagFieldEvaluate( const G4double y[], G4ThreeVector& Bfield );
      // Evaluate the field at a certain point.

  protected:  // without description

    // void MagFieldEvaluate( const G4double y[], G4double B[] )   
    //  { GetEquationOfMotion()->  GetFieldValue(y, B); }

  private:

    G4MagHelicalStepper(const G4MagHelicalStepper&);
    G4MagHelicalStepper& operator=(const G4MagHelicalStepper&);
      // Private copy constructor and assignment operator.

    static const G4double fUnitConstant;   //  As in G4Mag_EqRhs.hh/cc where it is not used.
  private:
  
    G4ThreeVector yInitial, yMidPoint, yFinal;
      // Data stored in order to find the chord.

    G4Mag_EqRhs*  fPtrMagEqOfMot;
};

#include  "G4MagHelicalStepper.icc"

#endif  /* G4MagHelicalStepper_hh */
