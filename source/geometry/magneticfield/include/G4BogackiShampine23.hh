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
// G4BogackiShampine23
//
// Class description:
//
// Bogacki-Shampine - 4 - 3(2) non-FSAL implementation 
//
// An implementation of the embedded RK method from the paper 
//  [1] P. Bogacki and L. F. Shampine,
//     "A 3(2) pair of Runge - Kutta formulas"
//     Appl. Math. Lett., vol. 2, no. 4, pp. 321-325, Jan. 1989.
//
// This version does not utilise the FSAL property of the method,
// which would allow the reuse of the last derivative in the next step.
// (Alternative FSAL implementation created with revised interface)

// Author: Somnath Banerjee (CERN, Google Summer of Code 2015), 20.05.2015
// Supervision: John Apostolakis (CERN)
// --------------------------------------------------------------------
#ifndef G4BOGACKI_SHAMPINE23_HH
#define G4BOGACKI_SHAMPINE23_HH

#include "G4MagIntegratorStepper.hh"
#include "G4FieldTrack.hh"

/**
 * @brief G4BogackiShampine23 is an integrator of particle's equation of
 * motion based on the Bogacki-Shampine non-FSAL implementation.
 */

class G4BogackiShampine23 : public G4MagIntegratorStepper
{
  public:

    /**
     * Constructor for G4BogackiShampine23.
     *  @param[in] EqRhs Pointer to the provided equation of motion.
     *  @param[in] numberOfVariables The number of integration variables.
     */
    G4BogackiShampine23(G4EquationOfMotion* EqRhs,
                        G4int numberOfVariables = 6);

    /**
     * Default Destructor.
     */
    ~G4BogackiShampine23() override = default;

    /**
     * Copy constructor and assignment operator not allowed.
     */
    G4BogackiShampine23(const G4BogackiShampine23&) = delete;
    G4BogackiShampine23& operator = (const G4BogackiShampine23&) = delete;

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
    void Stepper(const G4double yInput[],
                 const G4double dydx[],
                       G4double hstep,
                       G4double yOutput[],
                       G4double yError[]) override;

    /**
     * Same as the Stepper() function above, with dydx also in ouput.
     *  @param[in] yInput Starting values array of integration variables.
     *  @param[in] dydx Derivatives array.
     *  @param[in] hstep The given step size.
     *  @param[out] yOutput Integration output.
     *  @param[out] yError The estimated error.
     *  @param[out] dydxOutput dydx in output.
     */
    void Stepper(const G4double yInput[],
                 const G4double dydx[],
                       G4double hstep,
                       G4double yOutput[],
                       G4double yError[],
                       G4double dydxOutput[]);

    /**
     * Returns the distance from chord line.
     */
    G4double DistChord() const override;

    /**
     * Returns the order, 3, of integration.
     */
    inline G4int IntegratorOrder() const override { return 3; }

    /**
     * Returns the stepper type-ID, "kBogackiShampine23".
     */
    inline G4StepperType StepperType() const override { return kBogackiShampine23; }

  private:

    /**
     * Utility method used in Stepper() for computing the actual step.
     */
    void makeStep(const G4double yInput[],
                  const G4double dydx[],
                  const G4double hstep,
                        G4double yOutput[],
                        G4double* dydxOutput = nullptr,
                        G4double* yError = nullptr) const;

  private:

    G4double fyIn[G4FieldTrack::ncompSVEC],
             fdydx[G4FieldTrack::ncompSVEC],
             fyOut[G4FieldTrack::ncompSVEC],
             fdydxOut[G4FieldTrack::ncompSVEC];
    G4double fhstep = -1.0;
};

#endif
