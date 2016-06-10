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
// $Id: G4MagErrorStepper.hh 66356 2012-12-18 09:02:32Z gcosmo $
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
// --------------------------------------------------------------------

#ifndef G4MAGERRORSTEPPER_HH
#define G4MAGERRORSTEPPER_HH

#include "G4Types.hh"
#include "G4MagIntegratorStepper.hh"
#include "G4Mag_EqRhs.hh"
#include "G4ThreeVector.hh"

class G4MagErrorStepper : public G4MagIntegratorStepper
{
  public:  // with description

    G4MagErrorStepper(G4EquationOfMotion *EqRhs, G4int numberOfVariables, G4int numStateVariables=12);
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

    // STATE
    G4ThreeVector fInitialPoint, fMidPoint, fFinalPoint;
      // Data stored in order to find the chord
 
    // Dependent Objects, owned --- part of the STATE 
    G4double *yInitial, *yMiddle, *dydxMid, *yOneStep;
      // The following arrays are used only for temporary storage
      // they are allocated at the class level only for efficiency -
      // so that calls to new and delete are not made in Stepper().
};

#include  "G4MagErrorStepper.icc"

#endif  /* G4MAGERRORSTEPPER_HH */
