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

// Authors: Lucio Santi, Rodrigo Castro - 2018-2021.
// --------------------------------------------------------------------
#ifndef QSS_Stepper_HH
#define QSS_Stepper_HH

#include "G4FieldTrack.hh"
#include "G4FieldUtils.hh"
#include "G4LineSection.hh"
#include "G4MagIntegratorStepper.hh"
#include "G4QSS2.hh"
#include "G4QSS3.hh"
#include "G4QSSDriver.hh"
#include "G4QSSMessenger.hh"
#include "G4VIntegrationDriver.hh"
#include "G4qss_misc.hh"

#include <cmath>
#include <cassert>

// Maximum allowed number of QSS substeps per integration step
#define QSS_MAX_SUBSTEPS 1000

template <class QSS>
class G4QSStepper : public G4MagIntegratorStepper
{
  public:

    G4QSStepper(G4EquationOfMotion* EqRhs,
                G4int numberOfVariables = 6,
                G4bool primary = true);
   ~G4QSStepper() override;

    void Stepper(const G4double y[],
                 const G4double dydx[],
                       G4double h,
                       G4double yout[],
                       G4double yerr[]) override;

    void Stepper(const G4double yInput[],
                 const G4double dydx[],
                       G4double hstep,
                       G4double yOutput[],
                       G4double yError[],
                       G4double dydxOutput[]);

    // For calculating the output at the tau fraction of Step
    //
    inline void SetupInterpolation() {}
    inline void Interpolate(G4double tau, G4double yOut[]);

    G4double DistChord() const override;

    G4int IntegratorOrder() const override { return method->order(); }

    void reset(const G4FieldTrack* track);

    void SetPrecision(G4double dq_rel, G4double dq_min);
      // precision parameters for QSS method

    static G4QSStepper<G4QSS2>* build_QSS2(G4EquationOfMotion* EqRhs,
                                           G4int numberOfVariables = 6,
                                           G4bool primary = true);

    static G4QSStepper<G4QSS3>* build_QSS3(G4EquationOfMotion* EqRhs,
                                           G4int numberOfVariables = 6,
                                           G4bool primary = true);

    inline G4EquationOfMotion* GetSpecificEquation() { return GetEquationOfMotion(); }

    inline const field_utils::State& GetYOut() const { return fyOut; }

    inline G4double GetLastStepLength() { return fLastStepLength; }

  private:

    G4QSStepper(QSS* method,
                G4EquationOfMotion* EqRhs,
                G4int numberOfVariables = 6,
                G4bool primary = true);

    void initialize_data_structs();
    static QSS_simulator build_simulator();

    inline void update_field();
    inline void save_substep(G4double time, G4double length);

    inline void realloc_substeps();
    inline void get_state_from_poly(G4double* x, G4double* tx,
                                    G4double time, G4double* state);

    inline void recompute_derivatives(int index);
    inline void update_time();

    inline G4double get_coeff() { return fCoeff_local; }

    inline void set_coeff(G4double coeff) { fCoeff_local = coeff; }

    inline void set_charge(G4double q)
    {
      f_charge_c2 = q * cLight_local * cLight_local;  // 89875.5178737;
    }

    inline G4double get_qc2() { return f_charge_c2; }

    inline void set_mg() { fMassGamma = f_mass * fGamma2; }

    inline void set_gamma2(G4double gamma2) { fGamma2 = gamma2; }
    inline void set_velocity(G4double v) { fVelocity = v; }

    inline void velocity_to_momentum(G4double* state);

    inline void set_gamma(G4double p_sq)
    {
      set_gamma2(std::sqrt(p_sq / (f_mass * f_mass) + 1));
      set_mg();
      set_coeff(get_qc2() / fMassGamma);
    }

  private:

    QSS_simulator simulator;
    QSS* method;

    // State
    //
    G4double fLastStepLength;
    field_utils::State fyIn, fyOut;

    G4double f_mass;
    static constexpr G4double cLight_local = 299.792458;  // should use CLHEP
    G4double f_charge_c2;
    G4double fMassGamma;
    G4double fGamma2;
    G4double fCoeff_local;
    G4double fVelocity;
};

using G4QSStepper_QSS2 = G4QSStepper<G4QSS2>;
using G4QSStepper_QSS3 = G4QSStepper<G4QSS3>;

template <class QSS>
inline G4QSStepper<QSS>::G4QSStepper(QSS* qss, G4EquationOfMotion* EqRhs,
                              G4int noIntegrationVariables, G4bool)
  : G4MagIntegratorStepper(EqRhs, noIntegrationVariables),
    simulator(qss->getSimulator()),
    method(qss)
{
  SetIsQSS(true);  //  Replaces virtual method IsQSS
  fLastStepLength = -1.0;

  f_mass = 0;
  f_charge_c2 = 0;
  fMassGamma = 0;
  fGamma2 = 0;
  fCoeff_local = 0;
  fVelocity = 0;

  this->initialize_data_structs();
  this->SetPrecision(1e-4, 1e-7);  // Default values
}

template <class QSS>
inline G4QSStepper<QSS>::~G4QSStepper()
{
  for (auto & i : simulator->SD) { free(i); }

  free(SUBSTEPS(this->simulator));
  free(this->simulator);
}

template <class QSS>
inline void G4QSStepper<QSS>::Stepper(const G4double yInput[],
                                      const G4double dydx[],
                                            G4double hstep,
                                            G4double yOutput[],
                                            G4double yError[],
                                            G4double /*dydxOutput*/[])
{
  Stepper(yInput, dydx, hstep, yOutput, yError);
}

template <class QSS>
inline void G4QSStepper<QSS>::update_time()
{
  auto* const sim = this->simulator;

  sim->time = sim->nextStateTime[0];
  sim->minIndex = 0;

  if (sim->nextStateTime[1] < sim->time) {
    sim->time = sim->nextStateTime[1];
    sim->minIndex = 1;
  }
  if (sim->nextStateTime[2] < sim->time) {
    sim->time = sim->nextStateTime[2];
    sim->minIndex = 2;
  }
  if (sim->nextStateTime[3] < sim->time) {
    sim->time = sim->nextStateTime[3];
    sim->minIndex = 3;
  }
  if (sim->nextStateTime[4] < sim->time) {
    sim->time = sim->nextStateTime[4];
    sim->minIndex = 4;
  }
  if (sim->nextStateTime[5] < sim->time) {
    sim->time = sim->nextStateTime[5];
    sim->minIndex = 5;
  }
}

template <class QSS>
inline void G4QSStepper<QSS>::Stepper(const G4double yInput[],
                                      const G4double /*DyDx*/[],
                                            G4double max_length,
                                            G4double yOut[],
                                            G4double[] /*yErr[]*/)
{
  G4double elapsed;
  G4double t, prev_time = 0;
  G4double length = 0.;
  G4int index;

  const G4int coeffs = method->order() + 1;
  G4double* tq = simulator->tq;
  G4double* tx = simulator->tx;
  G4double* dQRel = simulator->dQRel;
  G4double* dQMin = simulator->dQMin;
  G4double* lqu = simulator->lqu;
  G4double* x = simulator->x;
  G4int** SD = simulator->SD;
  G4int cf0, infCf0;

  CUR_SUBSTEP(simulator) = 0;

  this->save_substep(0, length);

  this->update_time();
  t = simulator->time;
  index = simulator->minIndex;

  while (length < max_length && t < Qss_misc::INF && CUR_SUBSTEP(simulator) < QSS_MAX_SUBSTEPS) {
    cf0 = index * coeffs;
    elapsed = t - tx[index];
    method->advance_time_x(cf0, elapsed);
    tx[index] = t;
    lqu[index] = dQRel[index] * std::fabs(x[cf0]);
    if (lqu[index] < dQMin[index]) {
      lqu[index] = dQMin[index];
    }
    method->update_quantized_state(index);
    tq[index] = t;
    method->next_time(index, t);
    for (G4int i = 0; i < 3; i++) {
      G4int j = SD[index][i];
      elapsed = t - tx[j];
      infCf0 = j * coeffs;
      if (elapsed > 0) {
        x[infCf0] = method->evaluate_x_poly(infCf0, elapsed, x);
        tx[j] = t;
      }
    }

    this->update_field();
    this->recompute_derivatives(index);
    method->recompute_next_times(SD[index], t);

    if (t > prev_time) {
      length += fVelocity * (t - prev_time);
      if (length <= max_length) { this->save_substep(t, length); }
      else { break; }
    }

    this->update_time();
    prev_time = t;
    t = simulator->time;
    index = simulator->minIndex;
  }

  if(CUR_SUBSTEP(simulator) >= QSS_MAX_SUBSTEPS) {
    max_length = length;
  }

  auto* const substep = &LAST_SUBSTEP_STRUCT(simulator);
  t = substep->start_time + (max_length - substep->len) / fVelocity;

  this->get_state_from_poly(substep->x, substep->tx, t, yOut);

  velocity_to_momentum(yOut);

  const G4int numberOfVariables = GetNumberOfVariables();
  for (G4int i = 0; i < numberOfVariables; ++i) {
    // Store Input and Final values, for possible use in calculating chord
    fyIn[i] = yInput[i];
    fyOut[i] = yOut[i];
  }

  fLastStepLength = max_length;
}

template<class QSS>
inline G4double G4QSStepper<QSS>::DistChord() const
{
  G4double yMid[6];
  const_cast<G4QSStepper<QSS>*>(this)->Interpolate(0.5, yMid);

  const G4ThreeVector begin = makeVector(fyIn, field_utils::Value3D::Position);
  const G4ThreeVector end = makeVector(fyOut, field_utils::Value3D::Position);
  const G4ThreeVector mid = makeVector(yMid, field_utils::Value3D::Position);

  return G4LineSection::Distline(mid, begin, end);
}

template <class QSS>
inline void G4QSStepper<QSS>::Interpolate(G4double tau, G4double yOut[])
{
  G4double length = tau * fLastStepLength;
  G4int idx = 0, j = LAST_SUBSTEP(simulator);
  G4double end_time;

  if (j >= 15) {
    G4int i = 0, k = j;
    idx = j >> 1;
    while (idx < k && i < j - 1) {
      if (length < SUBSTEP_LEN(simulator, idx)) {
        j = idx;
      } else if (length >= SUBSTEP_LEN(simulator, idx + 1)) {
        i = idx;
      } else {
        break;
      }

      idx = (i + j) >> 1;
    }
  }
  else {
    for (; idx < j && length >= SUBSTEP_LEN(simulator, idx + 1); idx++) {;}
  }

  auto* const substep = &SUBSTEP_STRUCT(simulator, idx);
  end_time = substep->start_time + (length - substep->len) / fVelocity;

  this->get_state_from_poly(substep->x, substep->tx, end_time, yOut);

  velocity_to_momentum(yOut);
}

template <class QSS>
inline void G4QSStepper<QSS>::reset(const G4FieldTrack* track)
{
  using Qss_misc::PXidx;
  using Qss_misc::PYidx;
  using Qss_misc::PZidx;
  using Qss_misc::VXidx;
  using Qss_misc::VYidx;
  using Qss_misc::VZidx;

  G4ThreeVector pos = track->GetPosition();
  G4ThreeVector momentum = track->GetMomentum();

  f_mass = track->GetRestMass();
  set_charge(track->GetCharge());
  set_gamma(momentum.mag2());
  G4double c_mg = cLight_local / fMassGamma;
  set_velocity(momentum.mag() * c_mg);

  method->reset_state(PXidx, pos.getX());
  method->reset_state(PYidx, pos.getY());
  method->reset_state(PZidx, pos.getZ());

  method->reset_state(VXidx, momentum.getX() * c_mg);
  method->reset_state(VYidx, momentum.getY() * c_mg);
  method->reset_state(VZidx, momentum.getZ() * c_mg);

  this->update_field();
  method->full_definition(get_coeff());

  method->recompute_all_state_times(0);

  simulator->time = 0;
}

template <class QSS>
inline void G4QSStepper<QSS>::SetPrecision(G4double dq_rel, G4double dq_min)
{
  G4double* dQMin = simulator->dQMin;
  G4double* dQRel = simulator->dQRel;
  G4int n_vars = simulator->states;

  if (dq_min <= 0) { dq_min = dq_rel * 1e-3; }

  for (G4int i = 0; i < n_vars; ++i) {
    dQRel[i] = dq_rel;
    dQMin[i] = dq_min;
  }
}

template <class QSS>
inline void G4QSStepper<QSS>::initialize_data_structs()
{
  auto sim = this->simulator;
  auto  states = (G4int*)calloc(Qss_misc::VAR_IDX_END, sizeof(G4int));

  sim->states = Qss_misc::VAR_IDX_END;
  sim->it = 0.;

  for (unsigned int i = 0; i < Qss_misc::VAR_IDX_END; i++) {
    sim->SD[i] = (G4int*)malloc(3 * sizeof(G4int));
  }

  sim->SD[0][states[0]++] = 3;
  sim->SD[0][states[0]++] = 4;
  sim->SD[0][states[0]++] = 5;

  sim->SD[1][states[1]++] = 3;
  sim->SD[1][states[1]++] = 4;
  sim->SD[1][states[1]++] = 5;

  sim->SD[2][states[2]++] = 3;
  sim->SD[2][states[2]++] = 4;
  sim->SD[2][states[2]++] = 5;

  sim->SD[3][states[3]++] = 0;
  sim->SD[3][states[3]++] = 4;
  sim->SD[3][states[3]++] = 5;

  sim->SD[4][states[4]++] = 1;
  sim->SD[4][states[4]++] = 3;
  sim->SD[4][states[4]++] = 5;

  sim->SD[5][states[5]++] = 2;
  sim->SD[5][states[5]++] = 3;
  sim->SD[5][states[5]++] = 4;

  free(states);
}

template <class QSS>
inline QSS_simulator G4QSStepper<QSS>::build_simulator()
{
  QSS_simulator sim = (QSS_simulator)malloc(sizeof(*sim));
  MAX_SUBSTEP(sim) = Qss_misc::MIN_SUBSTEPS;
  SUBSTEPS(sim) = (QSSSubstep)malloc(Qss_misc::MIN_SUBSTEPS * sizeof(*SUBSTEPS(sim)));
  return sim;
}

template <class QSS>
inline void G4QSStepper<QSS>::recompute_derivatives(G4int index)
{
  const G4int coeffs = method->order() + 1;
  G4double e;
  G4int idx = 0;

  e = simulator->time - simulator->tq[0];
  if (likely(e > 0)) { method->advance_time_q(idx, e); }
  simulator->tq[0] = simulator->time;

  idx += coeffs;
  e = simulator->time - simulator->tq[1];
  if (likely(e > 0)) { method->advance_time_q(idx, e); }
  simulator->tq[1] = simulator->time;

  idx += coeffs;
  e = simulator->time - simulator->tq[2];
  if (likely(e > 0)) { method->advance_time_q(idx, e); }
  simulator->tq[2] = simulator->time;

  idx += coeffs;
  e = simulator->time - simulator->tq[3];
  if (likely(e > 0)) { method->advance_time_q(idx, e); }
  simulator->tq[3] = simulator->time;

  idx += coeffs;
  e = simulator->time - simulator->tq[4];
  if (likely(e > 0)) { method->advance_time_q(idx, e); }
  simulator->tq[4] = simulator->time;

  idx += coeffs;
  e = simulator->time - simulator->tq[5];
  if (likely(e > 0)) { method->advance_time_q(idx, e); }
  simulator->tq[5] = simulator->time;

  method->dependencies(index, get_coeff());
}

template <class QSS>
inline void G4QSStepper<QSS>::update_field()
{
  using Qss_misc::PXidx;
  using Qss_misc::PYidx;
  using Qss_misc::PZidx;

  const G4int order1 = method->order() + 1;
  G4double* const _field = simulator->alg;
  G4double* const _point = _field + order1;

  _point[PXidx] = simulator->x[PXidx];
  _point[PYidx] = simulator->x[PYidx * order1];
  _point[PZidx] = simulator->x[PZidx * order1];

  this->GetEquationOfMotion()->GetFieldValue(_point, _field);
}

template <class QSS>
inline void G4QSStepper<QSS>::save_substep(G4double time, G4double length)
{
  memcpy(CUR_SUBSTEP_X(simulator), simulator->x,
    (Qss_misc::VAR_IDX_END * (Qss_misc::MAX_QSS_STEPPER_ORDER + 2)) * sizeof(G4double));

  CUR_SUBSTEP_START(simulator) = time;
  CUR_SUBSTEP_LEN(simulator) = length;
  CUR_SUBSTEP(simulator)++;

  if (unlikely(CUR_SUBSTEP(simulator) == MAX_SUBSTEP(simulator))) {
    this->realloc_substeps();
  }
}

template <class QSS>
inline void G4QSStepper<QSS>::realloc_substeps()
{
  const G4int prev_index = MAX_SUBSTEP(simulator), new_index = 2 * prev_index;

  MAX_SUBSTEP(simulator) = new_index;
  SUBSTEPS(simulator) =
    (QSSSubstep)realloc(SUBSTEPS(simulator), new_index * sizeof(*SUBSTEPS(simulator)));
}

template <class QSS>
inline void G4QSStepper<QSS>::get_state_from_poly(
  G4double* x, G4double* tx, G4double time, G4double* state)
{
  unsigned int coeff_index = 0, i;
  const unsigned int x_order = method->order(), x_order1 = x_order + 1;

  for (i = 0; i < Qss_misc::VAR_IDX_END; ++i) {
    assert(tx[i] <= time);
    state[i] = method->evaluate_x_poly(coeff_index, time - tx[i], x);
    coeff_index += x_order1;
  }
}

template <class QSS>
inline void G4QSStepper<QSS>::velocity_to_momentum(G4double* state)
{
  using Qss_misc::VXidx;
  using Qss_misc::VYidx;
  using Qss_misc::VZidx;
  G4double coeff = fMassGamma / cLight_local;

  state[VXidx] *= coeff;
  state[VYidx] *= coeff;
  state[VZidx] *= coeff;
}

#endif
