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
// G4RK547FEq2
//
// Class description:
//
// An implementation of the 7 stage embedded Runge-Kutta 4,5 pair (RK547FEq2)
// from the paper:
//   D. J. Higham and G. Hall,
//   "Embedded Runge-Kutta formulae with stable equilibrium states",
//   J. Comput. Appl. Math., vol. 29, no. 1, pp. 25-33, 1990.

// Author: Dmitry Sorokin (CERN, Google Summer of Code 2017), 02.11.2017
// Supervision: John Apostolakis (CERN)
// --------------------------------------------------------------------
#ifndef G4RK547FEQ2_HH
#define G4RK547FEQ2_HH

#include "G4MagIntegratorStepper.hh"
#include "G4FieldTrack.hh"

/**
 * @brief G4RK547FEq2 is an implementation of the 7 stage embedded
 * Runge-Kutta 4,5 pair.
 */

class G4RK547FEq2 : public G4MagIntegratorStepper
{
  public:

    /**
     * Constructor for G4RK547FEq2.
     *  @param[in] EqRhs Pointer to the provided equation of motion.
     *  @param[in] integrationVariables The number of integration variables.
     */
    G4RK547FEq2(G4EquationOfMotion* EqRhs,
                G4int integrationVariables = 6);

    /**
     * Default Destructor.
     */
    ~G4RK547FEq2() override = default;

    /**
     * Copy constructor and assignment operator not allowed.
     */
    G4RK547FEq2 (const G4RK547FEq2&) = delete;
    G4RK547FEq2& operator = (const G4RK547FEq2&) = delete;

    /**
     * The stepper for the Runge Kutta integration.
     * The stepsize is fixed, with the step size given by 'hstep'.
     * Integrates ODE starting values yInput[0 to 6].
     * Outputs yOutput[] and its estimated error yError[].
     *  @param[in] yInput Starting values array of integration variables.
     *  @param[in] dydx Derivatives array.
     *  @param[in] hstep The given step size.
     *  @param[out] yOutput Integration output.
     *  @param[out] yError The estimated error.
     */
    void Stepper( const G4double yInput[],
                  const G4double dydx[],
                        G4double hstep,
                        G4double yOutput[],
                        G4double yError[] ) override;

    /**
     * Same as the Stepper() function above, with dydx also in ouput.
     *  @param[in] yInput Starting values array of integration variables.
     *  @param[in] dydx Derivatives array.
     *  @param[in] hstep The given step size.
     *  @param[out] yOutput Integration output.
     *  @param[out] yError The estimated error.
     *  @param[out] dydxOutput dydx in output.
     */
    void Stepper( const G4double yInput[],
                  const G4double dydx[],
                        G4double hstep,
                        G4double yOutput[],
                        G4double yError[],
                        G4double dydxOutput[] );

    /**
     * Returns the distance from chord line.
     */
    G4double DistChord() const override;

    /**
     * Returns the order, 4, of integration.
     */
    G4int IntegratorOrder() const override { return 4; }

    /**
     * Returns the stepper type-ID, "kRK547FEq2".
     */
    G4StepperType StepperType() const override { return kRK547FEq2; }

  private:

    /**
     * Utility method used in Stepper() for computing the actual step.
     */
    void makeStep( const G4double yInput[],
                   const G4double dydx[],
                   const G4double hstep,
                         G4double yOutput[],
                         G4double* dydxOutput = nullptr,
                         G4double* yError = nullptr ) const;

  private:

    G4double fyIn[G4FieldTrack::ncompSVEC],
             fdydx[G4FieldTrack::ncompSVEC],
             fyOut[G4FieldTrack::ncompSVEC],
             fdydxOut[G4FieldTrack::ncompSVEC];

    G4double fhstep= -1.0;
};

#endif
