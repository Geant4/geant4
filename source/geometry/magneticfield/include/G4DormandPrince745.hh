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
// G4DormandPrince745
//
// Class desription:
//
// An implementation of the 5th order embedded RK method from the paper:
// J. R. Dormand and P. J. Prince, "A family of embedded Runge-Kutta formulae"
// Journal of computational and applied Math., vol.6, no.1, pp.19-26, 1980.
//
// DormandPrince7 - 5(4) embedded RK method

// Author: Somnath Banerjee (CERN, Google Summer of Code 2015), 25.05.2015
// Supervision: John Apostolakis (CERN)
// --------------------------------------------------------------------
#ifndef G4DORMAND_PRINCE_745_HH
#define G4DORMAND_PRINCE_745_HH

#include "G4MagIntegratorStepper.hh"
#include "G4FieldUtils.hh"

/**
 * @brief G4DormandPrince745 implements the 5th order embedded Runge-Kutta
 * method, non-FSAL definition of the stepper() method that evaluates one step
 * in field propagation.
 */

class G4DormandPrince745 : public G4MagIntegratorStepper
{
  public:

    /**
     * Constructor for G4DormandPrince745.
     *  @param[in] equation Pointer to the provided equation of motion.
     *  @param[in] numberOfVariables The number of integration variables.
     */
    G4DormandPrince745(G4EquationOfMotion* equation,
                       G4int numberOfVariables = 6);

    /**
     * Default Destructor.
     */
    ~G4DormandPrince745() override = default;

    /**
     * Copy constructor and assignment operator not allowed.
     */
    G4DormandPrince745(const G4DormandPrince745&) = delete;
    G4DormandPrince745& operator=(const G4DormandPrince745&) = delete;
 
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
     *  @param[out] dydxOutput dysx in output.
     */
    void Stepper(const G4double yInput[],
                 const G4double dydx[],
                       G4double hstep,
                       G4double yOutput[],
                       G4double yError[],
                       G4double dydxOutput[]);


    /**
     * Interface method for interpolation setup. Does nothing here.
     */
    inline void SetupInterpolation() {}

    /**
     * Calculates the output at the tau fraction of Step.
     * Lower (4th) order interpolant given by Dormand and Prince.
     */
    void Interpolate4thOrder(G4double yOut[], G4double tau) const;

    /**
     * Wrapper for Interpolate4thOrder() function above.
     */
    inline void Interpolate(G4double tau, G4double yOut[]) const
    {
      Interpolate4thOrder(yOut, tau);
    }

    /**
     * Sets up the extra stages for the 5th order interpolant.
     */
    void SetupInterpolation5thOrder();

    /**
     * Calculates the interpolated result 'yOut' with the coefficients.
     * Interpolant of 5th order given by Baker, Dormand, Gilmore and Prince.
     */
    void Interpolate5thOrder(G4double yOut[], G4double tau) const;

    /**
     * Returns the distance from chord line.
     */
    G4double DistChord() const override;

    /**
     * Returns the order, 4, of integration.
     */
    inline G4int IntegratorOrder() const override { return 4; }

    /**
     * Returns the stepper type-ID, "kDormandPrince745".
     */
    inline G4StepperType StepperType() const override { return kDormandPrince745; }

    /**
     * Methods to return the stepper name and description.
     */
    const G4String& StepperTypeName() const;
    const G4String& StepperDescription() const;
   
    /**
     * Returns the field state in output.
     */
    inline const field_utils::State& GetYOut() const { return fyOut; }

    /**
     * Returns a pointer to the equation of motion.
     */
    inline G4EquationOfMotion* GetSpecificEquation() { return GetEquationOfMotion(); }

  private:

    field_utils::State ak2, ak3, ak4, ak5, ak6, ak7, ak8, ak9;
    field_utils::State fyIn, fyOut, fdydxIn;

    G4double fLastStepLength = -1.0;
};

#endif
