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
// G4SimpleHeum
//
// Class description:
//
// Simple Heum stepper for magnetic field:
//        x_1 = x_0  +
//              h * 1/4 * dx(t0,x0)  +
//                  3/4 * dx(t0+2/3*h, x0+2/3*h*(dx(t0+h/3,x0+h/3*dx(t0,x0)))) 
//
// Third order solver.

// Author: W.Wander (MIT), 12.09.1997
// -------------------------------------------------------------------
#ifndef G4SIMPLEHEUM_HH
#define G4SIMPLEHEUM_HH

#include "G4MagErrorStepper.hh"

/**
 * @brief G4SimpleHeum implements a simple Heum stepper for magnetic field
 * with 3rd order solver.
 */

class G4SimpleHeum : public G4MagErrorStepper
{
  public:

    /**
     * Constructor for G4SimpleHeum.
     *  @param[in] EqRhs Pointer to the provided equation of motion.
     *  @param[in] num_variables The number of integration variables.
     */
    G4SimpleHeum(G4EquationOfMotion* EqRhs,
                 G4int num_variables = 6);

    /**
     * Destructor.
     */
    ~G4SimpleHeum() override;

    /**
     * Copy constructor and assignment operator not allowed.
     */
    G4SimpleHeum(const G4SimpleHeum&) = delete;
    G4SimpleHeum& operator=(const G4SimpleHeum&) = delete;

    /**
     * The stepper for the Runge Kutta integration, but performing a 'dump' step
     * without error calculation.
     *  @param[in] y Starting values array of integration variables.
     *  @param[in] dydx Derivatives array.
     *  @param[in] h The given step size.
     *  @param[out] yout Integration output.
     */
    void DumbStepper( const G4double y[],
                      const G4double dydx[],
                            G4double h,
                            G4double yout[] ) override;
  
    /**
     * Returns the order, 3, of integration.
     */
    inline G4int IntegratorOrder() const override { return 3; }

    /**
     * Returns the stepper type-ID, "kSimpleHeum".
     */
    inline G4StepperType StepperType() const override { return kSimpleHeum; }

  private:

    G4int fNumberOfVariables = 0;

    G4double* dydxTemp = nullptr;
    G4double* dydxTemp2 = nullptr;
    G4double* yTemp = nullptr;
    G4double* yTemp2 = nullptr;
};

#endif
