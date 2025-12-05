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
// G4BFieldIntegrationDriver
//
// Class description:
//
// Specialised integration driver for pure magnetic field.

// Author: Dmitry Sorokin (CERN, Google Summer of Code 2017), 12.09.2018
// Supervision: John Apostolakis (CERN)
// --------------------------------------------------------------------
#ifndef G4BFIELD_INTEGRATION_DRIVER_HH
#define G4BFIELD_INTEGRATION_DRIVER_HH

#include "G4VIntegrationDriver.hh"
#include "G4Mag_EqRhs.hh"

#include <memory>

/**
 * @brief G4BFieldIntegrationDriver is specialised integration driver
 * for pure magnetic field.
 */

class G4BFieldIntegrationDriver : public G4VIntegrationDriver
{
  public:

    /**
     * Constructor for the integrator driver.
     *  @param[in] smallStepDriver Pointer to driver for small steps.
     *  @param[in] largeStepDriver Pointer to driver for large steps.
     */
    G4BFieldIntegrationDriver(
        std::unique_ptr<G4VIntegrationDriver> smallStepDriver, 
        std::unique_ptr<G4VIntegrationDriver> largeStepDriver);

    /**
     * Default Destructor.
     */
    ~G4BFieldIntegrationDriver() override = default;

    /**
     * Copy constructor and assignment operator not allowed.
     */
    G4BFieldIntegrationDriver(const G4BFieldIntegrationDriver &) = delete;
    const G4BFieldIntegrationDriver& operator =(const G4BFieldIntegrationDriver &) = delete;

    /**
     * Computes the step to take, based on chord limits.
     *  @param[in,out] track The current track in field.
     *  @param[in] hstep Proposed step length.
     *  @param[in] eps Requested accuracy, y_err/hstep.
     *  @param[in] chordDistance Maximum sagitta distance.
     *  @returns The length of step taken.
     */
    G4double AdvanceChordLimited(G4FieldTrack& track,
                                 G4double hstep,
                                 G4double eps,
                                 G4double chordDistance) override;

    /**
     * Integrates ODE from current s (s=s0) to s=s0+h with accuracy eps.
     *  @param[in,out] track The current track in field.
     *  @param[in] hstep Proposed step length.
     *  @param[in] eps Requested accuracy, y_err/hstep.
     *  @param[in] hinitial Initial minimum integration step.
     *  @returns true if integration succeeds.
     */
    inline G4bool AccurateAdvance(G4FieldTrack& track,
                                  G4double hstep,
                                  G4double eps,
                                  G4double hinitial = 0) override;

    /**
     * Checks whether the driver implements re-integration.
     *  @returns true if driver *Recalculates* when AccurateAdvance() is called.
     */
    inline G4bool DoesReIntegrate() const override;
   
    /**
     * Setter and getter for the equation of motion.
     */
    void SetEquationOfMotion(G4EquationOfMotion* equation) override;
    inline G4EquationOfMotion* GetEquationOfMotion() override;

    /**
     * Returns a pointer to the integrator stepper.
     */
    inline G4MagIntegratorStepper* GetStepper() override;

    /**
     * Computes a step size for the next step, taking the last step's
     * normalised error 'errMaxNorm'.
     *  @param[in] errMaxNorm The normalised error on last step.
     *  @param[in] hstepCurrent The current proposed step.
     *  @returns The step size for the next step.
     */
    inline G4double ComputeNewStepSize(G4double errMaxNorm,
                                       G4double hstepCurrent) override;

    /**
     * Setter and getter for verbosity.
     */
    inline void SetVerboseLevel(G4int level) override;
    inline G4int GetVerboseLevel() const override;

    /**
     * Dispatch interface method for computing step.
     */
    inline void OnComputeStep(const G4FieldTrack* track) override;

    /**
     * Dispatch interface method for initialisation/reset of driver.
     */
    inline void OnStartTracking() override;

    /**
     * Writes out to stream the parameters/state of the driver.
     */
    inline void StreamInfo( std::ostream& os ) const override;
   
    /**
     * Prints out statistics of the integrator driver.
     */
    void PrintStatistics() const;

    /** [[deprecated("will be removed")]] */
    inline void GetDerivatives(const G4FieldTrack& track,
                               G4double dydx[]) const override;

    /** [[deprecated("will be removed")]] */
    inline void GetDerivatives(const G4FieldTrack& track,
                               G4double dydx[],
                               G4double field[]) const override;

    /** [[deprecated("use GetEquationOfMotion() instead of GetStepper()->GetEquationOfMotion()")]] */
    inline const G4MagIntegratorStepper* GetStepper() const override;

  private:

    /**
     * Given the field track, computes the radius of the curvature in field.
     */
    G4double CurvatureRadius(const G4FieldTrack& track) const;

    /**
     * Returns the value of the field in the 'Field' array, give the track.
     *  @param[in] track The current field track.
     *  @param[in,out] Field The array with field values.
     */
    void GetFieldValue(const G4FieldTrack& track, 
                             G4double Field[] ) const;
   
  private:

    std::unique_ptr<G4VIntegrationDriver> fSmallStepDriver;
    std::unique_ptr<G4VIntegrationDriver> fLargeStepDriver;
    G4VIntegrationDriver* fCurrDriver = nullptr;
    G4Mag_EqRhs* fEquation = nullptr;

    G4int fSmallDriverSteps = 0;
    G4int fLargeDriverSteps = 0;
};

#include "G4BFieldIntegrationDriver.icc"

#endif
