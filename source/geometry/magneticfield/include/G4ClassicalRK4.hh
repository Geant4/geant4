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
// G4ClassicalRK4
//
// Class description:
//
// Integrate the equations of the motion of a particle in a magnetic field
// using the classical 4th Runge-Kutta method.

// Authors: J.Apostolakis, V.Grichine (CERN), 30.01.1997
// -------------------------------------------------------------------
#ifndef G4CLASSICALRK4_HH
#define G4CLASSICALRK4_HH

#include "G4MagErrorStepper.hh"

/**
 * @brief G4ClassicalRK4 integrates the equations of the motion of a particle
 * in a magnetic field using the classical 4th Runge-Kutta method.
 */

class G4ClassicalRK4 : public G4MagErrorStepper 
{
  public:

    /**
     * Constructor for G4ClassicalRK4.
     *  @param[in] EquationMotion Pointer to the provided equation of motion.
     *  @param[in] numberOfVariables The number of integration variables.
     */
    G4ClassicalRK4(G4EquationOfMotion* EquationMotion,
                   G4int numberOfVariables = 6) ;

    /**
     * Destructor.
     */
    ~G4ClassicalRK4() override ;

    /**
     * Copy constructor and assignment operator not allowed.
     */
    G4ClassicalRK4(const G4ClassicalRK4&) = delete;
    G4ClassicalRK4& operator=(const G4ClassicalRK4&) = delete;

    // A stepper that does not know about errors.
    // It is used by the MagErrorStepper stepper.
   
    /**
     * Given values for the variables y[0,..,n-1] and their derivatives
     * dydx[0,...,n-1] known at x, uses the classical 4th Runge-Kutta
     * method to advance the solution over an interval h and returns the
     * incremented variables as yout[0,...,n-1]. The user supplies the
     * function RightHandSide(x,y,dydx), which returns derivatives dydx at x.
     * The source is routine rk4 from NRC p.712-713.
     *  @param[in] yIn Starting values array of integration variables.
     *  @param[in] dydx Derivatives array.
     *  @param[in] h The given step size.
     *  @param[out] yOut Integration output.
     */
    void DumbStepper( const G4double yIn[],
                      const G4double dydx[],
                            G4double h,
                            G4double yOut[] ) override ;

    /**
     * Returns the order, 4, of integration.
     */
    G4int IntegratorOrder() const override { return 4; }

    /**
     * Returns the stepper type-ID, "kClassicalRK4".
     */
    G4StepperType StepperType() const override { return kClassicalRK4; }

  private:

    G4double *dydxm, *dydxt, *yt; // scratch space - not state 
};

#endif
