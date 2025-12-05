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
// G4MagIntegratorStepper
//
// Class description:
//
// Abstract base class for integrator of particle's equation of motion,
// used in tracking in space dependent magnetic field.
//
// A Stepper must integrate over NumberOfVariables elements,
// and also copy (from input to output) any of NoStateVariables  
// not included in the NumberOfVariables.  
// 
// So it is expected that NoStateVariables >= NumberOfVariables

// Author: John Apostolakis (CERN), 15.01.1997
// --------------------------------------------------------------------
#ifndef G4MAGINTEGRATORSTEPPER_HH
#define G4MAGINTEGRATORSTEPPER_HH

#include "G4Types.hh"
#include "G4EquationOfMotion.hh"
#include "G4FieldParameters.hh"
#include "G4VIntegrationDriver.hh"
#include "G4IntegrationDriver.hh"

class G4VIntegrationDriver;

/**
 * @brief G4MagIntegratorStepper is an abstract base class for integrator
 * of particle's equation of motion, used in tracking in space dependent
 * magnetic field.
 */

class G4MagIntegratorStepper
{
  public:

    /**
     * Constructor for G4MagIntegratorStepper.
     *  @param[in] Equation Pointer to the provided equation of motion.
     *  @param[in] numIntegrationVariables The number of integration variables.
     *  @param[in] numStateVariables The number of state variables.
     *  @param[in] isFSAL Flag to indicate if it is an FSAL (First Same As Last)
     *             type driver.
     */
    G4MagIntegratorStepper(G4EquationOfMotion* Equation, 
                           G4int numIntegrationVariables,
                           G4int numStateVariables = 12,
                           G4bool isFSAL = false );

    /**
     * Default virtual Destructor.
     */
    virtual ~G4MagIntegratorStepper() = default;

    /**
     * Copy constructor and assignment operator not allowed.
     */
    G4MagIntegratorStepper(const G4MagIntegratorStepper&) = delete;
    G4MagIntegratorStepper& operator=(const G4MagIntegratorStepper&) = delete;

    /**
     * The stepper for the Runge Kutta integration.
     * The stepsize is fixed, with the step size given by 'h'.
     * Integrates ODE starting values y[0 to 6].
     * Outputs yout[] and its estimated error yerr[].
     *  @param[in] y Starting values array of integration variables.
     *  @param[in] dydx Derivatives array.
     *  @param[in] h The given step size.
     *  @param[out] yout Integration output.
     *  @param[out] yerr The estimated error.
     */
    virtual void Stepper( const G4double y[],
                          const G4double dydx[],
                                G4double h,
                                G4double yout[],
                                G4double yerr[] ) = 0;

    /**
     * Estimates the maximum distance of a chord from the true path
     * over the segment last integrated.
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
    inline void RightHandSide( const G4double y[], G4double dydx[] ) const;

    /**
     * Calculates 'dydx' and 'field' at point 'y'.
     */
    inline void RightHandSide( const G4double y[],
                                     G4double dydx[],
                                     G4double field[] ) const;

    /**
     * Returns the number of variables that the stepper will integrate over.
     */
    inline G4int GetNumberOfVariables() const;

    /**
     * Returns the number of variables of state variables (>= above, integration).
     */
    inline G4int GetNumberOfStateVariables() const;

    /**
     * Returns the order of the integrator, i.e. its error behaviour is of
     * the order O(h^order).
     */
    virtual G4int IntegratorOrder() const = 0;

    /**
     * Returns the stepper type ID ('kUserStepper').
     * This function should be overriden in derived classes.
     */
    virtual G4StepperType StepperType() const { return kUserStepper; }

    /**
     * Replacement method - using new data member.
     */
    inline G4int IntegrationOrder();
   
    /**
     * Methods returning the pointer to the associated equation of motion.
     * As some steppers (e.g. RKG3) require other methods of Eq_Rhs this
     * function allows for access to them.
     */
    inline G4EquationOfMotion* GetEquationOfMotion();
    inline const G4EquationOfMotion* GetEquationOfMotion() const;

    /**
     * Setter for the equation of motion.
     */
    inline void SetEquationOfMotion(G4EquationOfMotion* newEquation); 

    /**
     * Methods for counting/resetting the number of calls to RHS method(s).
     */
    inline unsigned long GetfNoRHSCalls();
    inline void ResetfNORHSCalls();

    /**
     * Returns true if the stepper is of FSAL (First Same As Last) type.
     */
    inline G4bool IsFSAL() const;

    /**
     * Returns true if the stepper is of QSS (Quantum State Simulation) type.
     */
    inline G4bool isQSS() const;
    inline void SetIsQSS(G4bool val);

  protected:

    /**
     * Setters for the integration order and FSAL type.
     */
    inline void SetIntegrationOrder(G4int order);
    inline void SetFSAL(G4bool flag = true);
   
  private:

    G4EquationOfMotion* fEquation_Rhs = nullptr;
    const G4int fNoIntegrationVariables = 0; // Variables in integration
    const G4int fNoStateVariables = 0;       // Number required for FieldTrack

    /** Counter for calls to RHS method. */
    mutable unsigned long fNoRHSCalls = 0UL;

    // Parameters of a RK method -- must be shared by all steppers of a type
    // -- Invariants for a class
 
    G4int fIntegrationOrder = -1;  // must be set by stepper !!!
     // All ClassicalRK4 steppers are 4th order
    G4bool fIsFSAL = false;
     // Depends on RK method & implementation
    G4bool fIsQSS = false;
};

#include  "G4MagIntegratorStepper.icc"

#endif
