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
// G4QSSDriver
//
// QSS Interpolator Driver

// Authors: Lucio Santi, Rodrigo Castro (Univ. Buenos Aires), 2018-2021
// --------------------------------------------------------------------
#ifndef G4QSSDriver_HH
#define G4QSSDriver_HH

#include "G4InterpolationDriver.hh"
#include "G4QSSMessenger.hh"

/**
 * @brief G4QSSDriver is a templated driver class defining the QSS
 * (Quantum State Simulation) Interpolator Driver.
 */

template <class T>
class G4QSSDriver : public G4InterpolationDriver<T, true>
{
  public:

    /**
     * Constructor for G4QSSDriver.
     *  @param[in] T Pointer to the stepper algorithm.
     */
    inline G4QSSDriver(T* stepper);

    /**
     * Copy constructor and assignment operator not allowed.
     */
    G4QSSDriver(const G4QSSDriver&) = delete;
    const G4QSSDriver& operator=(const G4QSSDriver&) = delete;

    /**
     * Dispatch interface method for initialisation/reset of driver.
     */
    void OnStartTracking() override;

    /**
     * Computes the step to take, based on chord limits.
     *  @param[in,out] track The current track in field.
     *  @param[in] hstep Proposed step length.
     *  @param[in] eps Requested accuracy, y_err/hstep.
     *  @param[in] chordDistance Maximum sagitta distance.
     *  @returns The length of step taken.
     */
    inline G4double AdvanceChordLimited(G4FieldTrack& track,
                                        G4double hstep,
                                        G4double eps,
                                        G4double chordDistance) override;

    /**
     * Dispatch interface method for computing step.
     */
    inline void OnComputeStep(const G4FieldTrack* track) override;

    /**
     * Setter for driver precision parameters.
     */
    inline void SetPrecision(G4double dq_rel, G4double dq_min);

    /**
     * Takes one Step that is as large as possible while satisfying the
     * accuracy criterion.
     *  @param[in] it Stepper iterator.
     *  @param[in,out] y The current track state, y.
     *  @param[in] dydx dydx array.
     *  @param[in,out] hstep Step to attempt.
     *  @param[in] epsStep The relative accuracy.
     *  @param[in] curveLength Step start, x.
     *  @param[in,out] track Pointer to the Field track. Not used.
     *  @returns The step achieved.
     */
    inline G4double OneGoodStep(typename G4InterpolationDriver<T, true>::StepperIterator it,
                                field_utils::State& y,
                                field_utils::State& dydx,
                                G4double& hstep,
                                G4double epsStep,
                                G4double curveLength,
                                G4FieldTrack* track) override;

  private:

    using Base = G4InterpolationDriver<T, true>;

    G4bool initializedOnFirstRun = false;
};

#include "G4QSSDriver.icc"

#endif
