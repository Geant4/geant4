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
// $Id: G4CashKarpRKF45.hh 97598 2016-06-06 07:19:46Z gcosmo $
//
//
// class G4CashKarpRKF45
//
// Class description:
//
// The Cash-Karp Runge-Kutta-Fehlberg 4/5 method is an embedded fourth
// order method (giving fifth-order accuracy) for the solution of an ODE.
// Two different fourth order estimates are calculated; their difference
// gives an error estimate. [ref. Numerical Recipes in C, 2nd Edition]
// It is used to integrate the equations of the motion of a particle 
// in a magnetic field.

// History:
// - Created. J.Apostolakis, V.Grichine - 30.1.97
// -------------------------------------------------------------------

#ifndef G4CashKARP_RKF45
#define G4CashKARP_RKF45

#include "G4MagIntegratorStepper.hh"

class G4CashKarpRKF45 : public G4MagIntegratorStepper
{

  public:  // with description

    G4CashKarpRKF45( G4EquationOfMotion *EqRhs,
                     G4int numberOfVariables = 6,
                     G4bool primary= true ) ;
   ~G4CashKarpRKF45() ;

    void Stepper( const G4double y[],
                  const G4double dydx[],
                        G4double h,
                        G4double yout[],
                        G4double yerr[] ) ;

  public:  // without description

    G4double  DistChord()   const; 
    G4int IntegratorOrder() const { return 4; }

  private:

    void StepWithEst( const G4double yIn[],
                      const G4double dydx[],
                            G4double Step,
                            G4double yOut[],
                            G4double& alpha2,
                            G4double& beta2,
                      const G4double B1[],
                            G4double B2[]    );  
      // No longer used. Obsolete.

    G4CashKarpRKF45(const G4CashKarpRKF45&);
    G4CashKarpRKF45& operator=(const G4CashKarpRKF45&);
      // Private copy constructor and assignment operator.

  private:

   G4double *ak2, *ak3, *ak4, *ak5, *ak6, *yTemp, *yIn; // *ak7
      // scratch space

    G4double fLastStepLength;
    G4double *fLastInitialVector, *fLastFinalVector,
             *fLastDyDx, *fMidVector, *fMidError;
      // for DistChord calculations

    G4CashKarpRKF45* fAuxStepper; 

};

#endif /* G4CashKARP_RKF45 */
