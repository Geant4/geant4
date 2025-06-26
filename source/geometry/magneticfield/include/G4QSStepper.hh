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
// Authors - version 1 : Lucio Santi, Rodrigo Castro (Univ. Buenos Aires) - 2018-2021
//         - version 2 : Mattias Portnoy (Univ. Buenos Aires) - 2024
// --------------------------------------------------------------------

#ifndef G4QSS_STEPPER_HH
#define G4QSS_STEPPER_HH 1

#include "G4FieldTrack.hh"
#include "G4MagIntegratorStepper.hh"
#include "G4QSSMessenger.hh"
#include "G4QSSubstepStruct.hh"

#include <cmath>

class G4QSStepper : public G4MagIntegratorStepper
{

  public:

    G4QSStepper( G4EquationOfMotion* equation,
                 G4int num_integration_vars,
                 G4int num_state_vars,
                 G4bool isFSAL,
                 G4int verbosity=0 );

    G4QSStepper(G4EquationOfMotion *EqRhs,
                G4int numberOfVariables = 6,
                G4bool primary = true);

    virtual ~G4QSStepper();

    inline constexpr G4double Cubic_Function(const QSStateVector* states,
                                             G4int index, G4double delta_t);
   
    inline constexpr G4double Parabolic_Function(const QSStateVector* states,
                                                 G4int index, G4double delta_t);

    inline constexpr G4double Linear_Function(const QSStateVector* states,
                                              G4int index, G4double delta_t);

    /* 0 means position type, 1 means velocity type. */
    inline constexpr int INDEX_TYPE(G4int i);

    inline void set_qss_order(G4int order);

    // auxiliary methods

    inline void momentum_to_velocity(const G4double* momentum, G4double* out);

    void set_relativistic_coeff(const G4double* momentum);

    inline void velocity_to_momentum(G4double *y);

    // Key methods

    void initialize(const G4double y[]);

    inline void compare_time_and_update(G4int index, G4int i);

    inline G4int get_next_sync_index();

    inline void update_field();

    inline G4double extrapolate_polynomial(QSStateVector* states,
                                    G4int index, G4double delta_t, G4int order);
    inline void extrapolate_all_states_to_t(Substep* substep,
                                            G4double t, G4double* yOut);

    /* Moves all the x states of variable index to the current time t. */
    inline void update_x(G4int index, G4double t);

    /* Moves all the q states of variable index to the current t. */
    inline void update_q(G4int index, G4double t);

    inline void update_x_position_derivates_using_q(G4int index);
    inline void update_x_velocity_derivates_using_q(G4int index);
    inline void update_x_derivates_using_q(G4int index);
    inline void update_sync_time_one_coefficient(G4int index);

    /* Updates when does the x,q distance goes beyond the quantum.
       Uses polynomial roots-finding formulas. */
    void update_sync_time(G4int index);

    /*  Key method called by driver. */
    void Stepper( const G4double y[],
                  const G4double /*dydx*/ [],
                  G4double h,
                  G4double yout[],
                  G4double /* yerr */ [] ) override;

    /* Obligatory G4InterpolationDriver methods. */
    inline G4int IntegratorOrder() const override;
    inline G4EquationOfMotion* GetSpecificEquation();
    inline const field_utils::State& GetYOut() const;

    void Interpolate(G4double tau,G4double yOut[]);

    inline G4double DistChord() const override;

    inline void Stepper(const G4double yInput[],
                        const G4double dydx[],
                        G4double hstep, G4double yOutput[], G4double yError[],
                        G4double /*dydxOutput*/ []);

    inline void SetupInterpolation();

    /* obligatory qss driver methods. */

    inline void reset(const G4FieldTrack* track);

    inline void SetPrecision(G4double dq_rel, G4double dq_min);

    inline G4double GetLastStepLength();

  private:

    // Constants

    static constexpr int DERIVATIVE_0 = 0;
    static constexpr int DERIVATIVE_1 = 1;
    static constexpr int DERIVATIVE_2 = 2;
    static constexpr int DERIVATIVE_3 = 3;
   
    static constexpr int VX = 3;
    static constexpr int VY = 4;
    static constexpr int VZ = 5;

    static constexpr int POSITION_IDX = 0;
    static constexpr int VELOCITY_IDX = 3;
    static constexpr int NUMBER_OF_VARIABLES_QSS = 6;

    static constexpr G4double INFTY = 1e+20;

    /* Used to check if field changed from last update field during substeps. */
    G4bool fField_changed = true;
    G4bool fTrack_changed = true;

    G4int  qss_order = 2;
   
    Substeps substeps;
    Substep current_substep;
    const G4FieldTrack* fCurrent_track = nullptr;
    QSStateVector dq_vector;

    // Invariants for this track -- during propagation
    //
    G4double fCharge;
    G4double fCharge_c2;
    G4double fRestMass;
    G4double fGamma;
    G4double fCoeff; // coeff;

    // Cached values -- for tiny speed up
    //
    G4double fMassOverC ; // was mass_times_gamma_over_speed_of_light;
    G4double fInv_mass_over_c;

    /* used by interpolation driver, need to copy state here
       when stepper finished. */
    G4double fYout[12];

    // QSS parameters separated into velocity and position
    //
    G4double dqrel[2] = {0.0,0.0};
    G4double dqmin[2] = {0.001,0.001};

    G4double fVelocity;
    G4double fFinal_t;
};

// ----------------------------------------------------------------------------
// Inline methods
// ----------------------------------------------------------------------------

#include "G4QSStepper.icc"

#endif
