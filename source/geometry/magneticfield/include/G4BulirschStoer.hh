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
// G4BulirschStoer
//
// Class description:
//
// The Bulirsch-Stoer is a controlled driver that adjusts both step size
// and order of the method. The algorithm uses the modified midpoint and
// a polynomial extrapolation computes the solution.

// Author: Dmitry Sorokin (CERN, Google Summer of Code 2016), 13.02.2018
// Supervision: John Apostolakis (CERN)
// --------------------------------------------------------------------
#ifndef G4BULIRSCH_STOER_HH
#define G4BULIRSCH_STOER_HH

#include "G4ModifiedMidpoint.hh"

#include "G4FieldTrack.hh"

/**
 * @brief G4BulirschStoer is a controlled driver that adjusts both step size
 * and order of the method. The algorithm uses the modified midpoint and
 * a polynomial extrapolation computes the solution.
 */

class G4BulirschStoer
{
  public:

    enum class step_result { success, fail };

    /**
     * Constructor for G4BulirschStoer.
     *  @param[in] equation Pointer to the provided equation of motion.
     *  @param[in] nvar The number of integration variables.
     *  @param[in] eps_rel Relative tolerance.
     *  @param[in] max_dt Maximum allowed time step.
     */
    G4BulirschStoer(G4EquationOfMotion* equation, G4int nvar,
                    G4double eps_rel, G4double max_dt = DBL_MAX);

    /**
     * Default Destructor.
     */
    ~G4BulirschStoer() = default;

    /**
     * Modifiers.
     */
    inline void set_max_dt(G4double max_dt);
    inline void set_max_relative_error(G4double eps_rel);

    /**
     * Stepper method.
     *  @param[in] in Initial position.
     *  @param[in] dxdt dxdt for mid-point calculation.
     *  @param[out] t The updated step.
     *  @param[out] out Updated position.
     *  @param[in,out] dt Step size.
     *  @returns success if step is not rejected.
     */
    step_result try_step(const G4double in[], const G4double dxdt[],
                         G4double& t, G4double out[], G4double& dt);

    /**
     * Resets the internal state of the stepper.
     */
    void reset();

    /**
     * Setter and getter for the equation of motion.
     */
    inline void SetEquationOfMotion(G4EquationOfMotion* equation);
    inline G4EquationOfMotion* GetEquationOfMotion() const;

    /**
     * Returns the number of integration variables.
     */
    inline G4int GetNumberOfVariables() const;

  private:

    /**
     * Polynomial extrapolation.
     */
    void extrapolate(std::size_t k, G4double xest[]);

    /**
     * Calculates the optimal step size for a given error and stage number.
     */
    G4double calc_h_opt(G4double h, G4double error, std::size_t k) const;

    /**
     * Calculates the optimal stage number.
     */
    G4bool set_k_opt(std::size_t k, G4double& dt);

    /**
     * Utilities.
     */
    G4bool in_convergence_window(G4int k) const;
    G4bool should_reject(G4double error, G4int k) const;

  private:

    /** Maximum number of stages. */
    const static G4int m_k_max = 8;

    /** Number of vars to be integrated. */
    G4int fnvar;

    /** Relative tolerance. */
    G4double m_eps_rel;

    /** Modified midpoint algorithm. */
    G4ModifiedMidpoint m_midpoint;

    /** Flags for step. */
    G4bool m_last_step_rejected{false};
    G4bool m_first{true};

    /** Last step size. */
    G4double m_dt_last{0.0};

    /** Max allowed time step. */
    G4double m_max_dt;

    /** Crude estimate of optimal order. */
    G4int m_current_k_opt;

    /** Error estimate. */
    G4double m_err[G4FieldTrack::ncompSVEC];

    /** Stores the successive interval counts. */
    G4int m_interval_sequence[m_k_max+1];

    /** Extrapolation coeffs (Neville's algorithm). */
    G4double m_coeff[m_k_max+1][m_k_max];

    /** Costs for interval count. */
    G4int m_cost[m_k_max+1];

    /** Sequence of states for extrapolation. */
    G4double m_table[m_k_max][G4FieldTrack::ncompSVEC];

    /** Optimal step size. */
    G4double h_opt[m_k_max+1];

    /** Work per unit step. */
    G4double work[m_k_max+1];
};

#include "G4BulirschStoer.icc"

#endif
