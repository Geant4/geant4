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
// G4VFSALIntegrationStepper
//
// Class description:
//
// Class similar to G4VMagIntegratorStepper, for steppers which
// estimate the value of the derivative at the projected endpoint
// of integration - at each successful step.
// This ability is known as 'First Same As Last' (FSAL). It
// reduces the number of required calls to the equation's 
// RightHandSide method, and, as such the number of calls to the 
// (potentially expensive) field evaluation methods.
//
// Based on G4VMagIntegratorStepper

// Author: Somnath Banerjee (CERN, Google Summer of Code 2015), 26.05.2016
// Supervision: John Apostolakis (CERN)
// --------------------------------------------------------------------
#ifndef G4VFSALINTEGRATOR_STEPPER_HH
#define G4VFSALINTEGRATOR_STEPPER_HH

#include "G4Types.hh"
#include "G4EquationOfMotion.hh"

/**
 * @brief G4VFSALIntegrationStepper is a class similar to
 * G4VMagIntegratorStepper, but for steppers which estimate the value of
 * the derivative at the projected endpoint of integration, at each successful
 * step. This ability is known as 'First Same As Last' (FSAL).
 * It reduces the number of required calls to the equation's RightHandSide
 * method, and, as such the number of calls to the (potentially expensive)
 * field evaluation methods.
 */

class G4VFSALIntegrationStepper
{
  public:

    /**
     * Constructor for G4VFSALIntegrationStepper.
     *  @param[in] Equation Pointer to the provided equation of motion.
     *  @param[in] numStateVariables The number of state variables.
     */
    G4VFSALIntegrationStepper (G4EquationOfMotion* Equation,
                               G4int numIntegrationVariables,
                               G4int numStateVariables = 12);

    /**
     * Default Destructor.
     */
    virtual ~G4VFSALIntegrationStepper() = default;

    /**
     * Copy constructor and assignment operator not allowed.
     */
    G4VFSALIntegrationStepper(const G4VFSALIntegrationStepper&) = delete;
    G4VFSALIntegrationStepper& operator=(const G4VFSALIntegrationStepper&) = delete;

    /**
     * The stepper for the Runge Kutta integration.
     * The stepsize is fixed, with the step size given by 'h'.
     * Integrates ODE starting values yInput[0 to 6].
     * Outputs yout[] and its estimated error yerr[].
     *  @param[in] y Starting values array of integration variables.
     *  @param[in] dydx Derivatives array.
     *  @param[in] h The given step size.
     *  @param[out] yout Integration output.
     *  @param[out] yerr The estimated error.
     *  @param[out] lastDydx Last derivative.
     */
    virtual void Stepper( const G4double y[],
                          const G4double dydx[],
                                G4double h,
                                G4double yout[],
                                G4double yerr[],
                                G4double lastDydx[]) = 0;

    /**
     * Returns an estimate of the maximum distance of a chord from the
     * true path over the segment last integrated.
     */
    virtual G4double DistChord() const = 0; 
    
    /**
     * Simple utility function to (re)normalise 'unit velocity' vector.
     */
    inline void NormaliseTangentVector( G4double vec[6] );

    /**
     * Simple utility function to (re)normalise 'unit spin' vector.
     */
    inline void NormalisePolarizationVector( G4double vec[12] );

    /**
     * Utility method to supply the standard Evaluation of the
     * Right Hand side of the associated equation.
     */
    void RightHandSide( const double y[], double dydx[] );

    /**
     * Returns the number of variables that the stepper will integrate over.
     */
    inline G4int GetNumberOfVariables() const;

    /**
     * Returns the number of variables of state variables (>= above, integration)
     */
    inline G4int GetNumberOfStateVariables() const;

    /**
     * Returns the order of the integrator, i.e. its error behaviour
     * is of the order O(h^order).
     */
    virtual G4int IntegratorOrder() const = 0;

    /**
     * Returns a pointer to the equation of motion.
     * As some steppers (e.g. RKG3) require other methods of Eq_Rhs,
     * this function allows for access to them.
     */
    inline G4EquationOfMotion* GetEquationOfMotion(); 

    /**
     * Setter for the equation of motion.
     */
    inline void SetEquationOfMotion(G4EquationOfMotion* newEquation); 

    /**
     * Methods for debug use.
     */
    inline G4int GetfNoRHSCalls() { return fNoRHSCalls; }
    void increasefNORHSCalls();
    inline void ResetfNORHSCalls() { fNoRHSCalls = 0; }

  private:

    G4EquationOfMotion* fEquation_Rhs = nullptr;

    /** Variables in integration. */
    const G4int fNoIntegrationVariables = 0;

    /** Number required for FieldTrack. */
    const G4int fNoStateVariables = 0;

    /** Used for debug. */
    G4int fNoRHSCalls = 0;
};

#include "G4VFSALIntegrationStepper.icc"

#endif
