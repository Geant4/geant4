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
// $Id: G4MagErrorStepper.hh,v 1.8 2001-07-11 09:59:08 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
//
// class G4MagErrorStepper
//
// Class description:
//
// Abstract base class for integrator of particle's equation of motion,
// used in tracking in space dependent magnetic field.

// History:
// 09.12.97  W.Wander <wwc@mit.edu>  Created G4MagErrorStepper
// 09.03.98  W.Wander <wwc@mit.edu>  Added AdvanceHelix functionality
// 09.11.98  J.Apostolakis           Moved AdvanceHelix to G4MagHelicalStepper

#ifndef G4MAGERRORSTEPPER_HH
#define G4MAGERRORSTEPPER_HH

#include "globals.hh"
#include "G4MagIntegratorStepper.hh"
#include "G4Mag_EqRhs.hh"
#include "G4ThreeVector.hh"

class G4MagErrorStepper : public G4MagIntegratorStepper
{
  public:  // with description

    G4MagErrorStepper(G4Mag_EqRhs *EqRhs,G4int numberOfVariables);
    virtual ~G4MagErrorStepper();
  
    void Stepper( const G4double y[],
		  const G4double dydx[],
		        G4double h,
		        G4double yout[],
		        G4double yerr[]  );
      // The stepper for the Runge Kutta integration. The stepsize 
      // is fixed, with the Step size given by h.
      // Integrates ODE starting values y[0 to 6].
      // Outputs yout[] and its estimated error yerr[].

    virtual  void DumbStepper( const G4double y[],
			       const G4double dydx[],
			             G4double h,
			             G4double yout[] ) = 0;
      // Performs a 'dump' Step without error calculation.

    G4double DistChord() const;

  private:

    G4MagErrorStepper(const G4MagErrorStepper&);
    G4MagErrorStepper& operator=(const G4MagErrorStepper&);
      // Private copy constructor and assignment operator.

  private:

    G4ThreeVector fInitialPoint, fMidPoint, fFinalPoint;
      // Data stored in order to find the chord
 
    // G4int theNumberOfVariables ; 

    G4double *yInitial, *yMiddle, *dydxMid, *yOneStep;
      // The following arrays are used only for temporary storage
      // they are allocated at the class level only for efficiency -
      // so that calls to new and delete are not made in Stepper().
};

#include  "G4MagErrorStepper.icc"

#endif  /* G4MAGERRORSTEPPER_HH */
