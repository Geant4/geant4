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
// G4QSStepper
//
// QSS Integrator Stepper
//
// Authors: version 1 - Lucio Santi, Rodrigo Castro (Univ. Buenos Aires), 2018-2021
//          version 2 - Mattias Portnoy (Univ. Buenos Aires), 2024
// --------------------------------------------------------------------
#ifndef G4QSS_STEPPER_HH
#define G4QSS_STEPPER_HH

#include "G4FieldTrack.hh"
#include "G4MagIntegratorStepper.hh"
#include "G4QSSubstepStruct.hh"

#include <cmath>
#include <CLHEP/Units/PhysicalConstants.h>

/**
 * @brief G4QSStepper is an integrator of particle's equation of
 * motion based on the QSS implementation.
 */

class G4QSStepper : public G4MagIntegratorStepper
{
  public:

    /**
     * Constructor for G4QSStepper.
     *  @param[in] equation Pointer to the provided equation of motion.
     *  @param[in] num_integration_vars The number of integration variables.
     *  @param[in] qssOrder The QSS order (2 or 3 expected; if <= 0 , use value
     *             from Messenger.
     */
    G4QSStepper( G4EquationOfMotion* equation,
                 G4int num_integration_vars = 6,  // always 6 -- ignore
                 G4int qssOrder= -1 ); 

    /**
     * Default Destructor. Freeing of memory is done in susbsteps destructor.
     */
    ~G4QSStepper() override = default;

    /**
     * Utility methods.
     */
    inline constexpr G4double Cubic_Function(const QSStateVector* states,
                                             G4int index, G4double delta_t);
    inline constexpr G4double Parabolic_Function(const QSStateVector* states,
                                                 G4int index, G4double delta_t);
    inline constexpr G4double Linear_Function(const QSStateVector* states,
                                              G4int index, G4double delta_t);

    /**
     * 0 means position type, 1 means velocity type.
     */
    inline constexpr int INDEX_TYPE(G4int i);

    /**
     * Auxiliary methods.
     */
    inline void momentum_to_velocity(const G4double* momentum, G4double* out);
    void set_relativistic_coeff(const G4double* momentum);
    inline void velocity_to_momentum(G4double *y);

    /**
     * Key methods.
     */
    void initialize(const G4double y[]);
    inline void compare_time_and_update(G4int& index, G4int i);
    inline G4int get_next_sync_index();
    inline void update_field();
    inline G4double extrapolate_polynomial(QSStateVector* states,
                                    G4int index, G4double delta_t, G4int order);
    inline void extrapolate_all_states_to_t(Substep* substep,
                                            G4double t, G4double* yOut);

    /**
     * Moves all the x states of variable index to the current time t.
     */
    inline void update_x(G4int index, G4double t);

    /**
     * Moves all the q states of variable index to the current t.
     */
    inline void update_q(G4int index, G4double t);

    /**
     * Update methods.
     */
    inline void update_x_position_derivates_using_q(G4int index);
    inline void update_x_velocity_derivates_using_q(G4int index);
    inline void update_x_derivates_using_q(G4int index);
    inline void update_sync_time_one_coefficient(G4int index);

    /*
     * Updates when does the x,q distance goes beyond the quantum.
     * Uses polynomial roots-finding formulas.
     */
    void update_sync_time(G4int index);

    /**
     * The stepper for the integration.
     * The stepsize is fixed, with the step size given by 'h'.
     * Integrates ODE starting values y[0 to 6]. Outputs yout[].
     *  @param[in] y Starting values array of integration variables.
     *  @param[in] dydx Derivatives array - Not used.
     *  @param[in] h The given step size.
     *  @param[out] yout Integration output.
     *  @param[out] yError The estimated error - Not used.
     */
    void Stepper( const G4double y[],
                  const G4double /*dydx*/ [],
                  G4double h,
                  G4double yout[],
                  G4double /* yerr */ [] ) override;

    /**
     * Returns the QSS order of integration.
     */
    inline G4int IntegratorOrder() const override;

    /**
     * Returns the stepper type-ID, "kQSStepper".
     */
    inline G4StepperType StepperType() const override { return kQSStepper; }

    /**
     * Returns a pointer to the equation of motion.
     */
    inline G4EquationOfMotion* GetSpecificEquation();

    /**
     * Returns current track state.
     */
    inline const field_utils::State& GetYOut() const;

    /**
     * Track interpolation.
     *  @param[in] tau Step start, x.
     *  @param[in,out] yOut The current track state, y.
     */
    void Interpolate(G4double tau, G4double yOut[]);

    /**
     * Returns the distance from chord line.
     */
    inline G4double DistChord() const override;

    /**
     * Wrapper for the Stepper() function above.
     */
    inline void Stepper(const G4double yInput[],
                        const G4double dydx[],
                        G4double hstep,
                        G4double yOutput[],
                        G4double yError[],
                        G4double /*dydxOutput*/ []);

    /**
     * Sets up interpolation. Does nothing.
     */
    inline void SetupInterpolation();

    /*
     * Obligatory qss driver methods.
     */
    inline void reset(const G4FieldTrack* track);
    inline void SetPrecision(G4double dq_rel, G4double dq_min);
    inline G4double GetLastStepLength();

    /*
     * Sets the mass at rest. Checking/ensuring that it is positive.
     */
    inline void setRestMass(G4double restMass); 

  private:

    // Constants

    static constexpr int DERIVATIVE_0{0};
    static constexpr int DERIVATIVE_1{1};
    static constexpr int DERIVATIVE_2{2};
    static constexpr int DERIVATIVE_3{3};
   
    static constexpr int VX{3};
    static constexpr int VY{4};
    static constexpr int VZ{5};

    static constexpr int POSITION_IDX{0};
    static constexpr int VELOCITY_IDX{3};

    static constexpr G4double INFTY{1e+20};

    /** Used to check if field changed from last update field during substeps. */
    G4bool fField_changed{true};
    G4bool fTrack_changed{true};

    const G4int qss_order{2};
   
    Substeps substeps;
    Substep current_substep;
    const G4FieldTrack* fCurrent_track{nullptr};
    QSStateVector dq_vector;

    /** Invariants for this track -- during propagation. */
    G4double fCharge{-1.0};
    G4double fCharge_c2;
    G4double fRestMass{CLHEP::electron_mass_c2};
    G4double fGamma{1.0};
    G4double fCoeff; // coeff;

    /** Cached values -- for tiny speed up. */
    G4double fMassOverC ; // was mass_times_gamma_over_speed_of_light;
    G4double fInv_mass_over_c;

    /** Used by interpolation driver, need to copy state here when stepper finished. */
    G4double fYout[12];

    /** QSS parameters separated into velocity and position. */
    G4double dqrel[2] = {0.0,0.0};
    G4double dqmin[2] = {0.001,0.001};

    G4double fVelocity{0.0};
    G4double fFinal_t{0.0};
};

// ----------------------------------------------------------------------------
// Inline methods
// ----------------------------------------------------------------------------

#include "G4QSStepper.icc"

#endif
