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

#include "G4QSStepper.hh"
#include "G4QSSMessenger.hh"
#include "G4AutoLock.hh"

#include <algorithm>

namespace
{
  // Mutex to lock updating the global ion map
  G4Mutex qssOrderCheckMutex = G4MUTEX_INITIALIZER;
}

// ----------------------------------------------------------------------------

G4QSStepper::G4QSStepper( G4EquationOfMotion* equation,
                          G4int   , // num_integration_vars is always 6 -- Needed for ctor in Integration Driver)
                          G4int  qssOrder ):
  G4MagIntegratorStepper(equation, NUMBER_OF_VARIABLES_QSS, NUMBER_OF_VARIABLES_QSS, false), 
                                // num_integration_vars,    num_state_vars,      isFSAL)
  qss_order( qssOrder <= 0 ? G4QSSMessenger::instance()->GetQssOrder() : qssOrder ) 
{
  using std::memset;
  SetIsQSS(true);

  // Note: current_substep is initialised by 'Substep' constructor

  if( qss_order < 2 || qss_order > 3 )
  {
     G4ExceptionDescription err_msg;     
     err_msg << "-G4QSStepper/c-tor: qss_order= " << qss_order << "  is not either 2 or 3 ";
     G4Exception("G4QSStepper::G4QSStepper", "GeomMag-0001", FatalException, err_msg );
  }

  auto messenger= G4QSSMessenger::instance();

  if( qss_order != messenger->GetQssOrder() ) 
  {
    // BARRIER
    G4AutoLock lock(qssOrderCheckMutex);

    if( qss_order != messenger->GetQssOrder() ) 
    {
      static std::atomic<G4int> change_of_order= 0;
      change_of_order++;

      if( change_of_order == 1 )
      {
        G4QSSMessenger::instance()->SetQssOrder(qss_order);
      }
      else
      {
        G4ExceptionDescription err_msg;
        err_msg << "G4QSSStepper: Trying to change order of QSS stepper(s) for " 
                << change_of_order << " th time -- Maximum allowed is 1 time. "
                << " Requesting order = " << qss_order 
                << " whereas current (static) value = " << messenger->GetQssOrder()
                << G4endl;
        G4Exception("G4QSStepper::G4QSStepper","G4QSS-0002",FatalException, err_msg);
      }
    }
  }

  // Place-holder values -- will be initialised for each particle track
  fCharge_c2 = fCharge * ( 1.0 / (CLHEP::c_light*CLHEP::c_light) ); 
  G4double momentum[3]= { 0.0, 0.0, 0.0 };
  set_relativistic_coeff(momentum);
  std::fill_n( fYout, 12, 0.0 );
}


// ----------------------------------------------------------------------------

void G4QSStepper::set_relativistic_coeff(const G4double* momentum) 
{
  G4double momentum2 = momentum[0]*momentum[0] + momentum[1]*momentum[1] + momentum[2]*momentum[2];
  fGamma = std::sqrt(momentum2/(fRestMass*fRestMass) + 1);
  G4double mass_times_gamma = fRestMass * fGamma;
  fMassOverC = mass_times_gamma * (1.0 / CLHEP::c_light);
  fInv_mass_over_c = CLHEP::c_light * (1.0 / mass_times_gamma);
  fCoeff = fCharge_c2 / mass_times_gamma;
}

// ----------------------------------------------------------------------------

void G4QSStepper::initialize(const G4double y[]) 
{
  using std::memset;

  substeps.reset();

  // Load values particle's invariants
  if (fCurrent_track != nullptr)
  {
    fCharge = fCurrent_track->GetCharge();
    fCharge_c2 = fCharge * CLHEP::c_squared;    /// Was 89875.5178737;
    setRestMass( fCurrent_track->GetRestMass() );
  }

  // y contains postion in first 3 index and momentum on the next 3
  set_relativistic_coeff(&y[3]);

  G4double velocity_vector[3];
  momentum_to_velocity(&y[3], velocity_vector);
  fVelocity = std::sqrt(velocity_vector[0]*velocity_vector[0] + velocity_vector[1]*velocity_vector[1] + velocity_vector[2]*velocity_vector[2] );

  std::copy_n( y,              3,  &current_substep.state_x[DERIVATIVE_0][POSITION_IDX] );
  std::copy_n( velocity_vector, 3, &current_substep.state_x[DERIVATIVE_0][VELOCITY_IDX] );
  // Copy into state_q
  std::copy_n( &current_substep.state_x[DERIVATIVE_0][0], 6, &current_substep.state_q[DERIVATIVE_0][0] );

  for (G4int i = 1; i < qss_order; ++i)
  {
    std::fill_n(current_substep.state_q[i], NUMBER_OF_VARIABLES_QSS, 0.0);
  }

  std::fill_n(current_substep.state_tx, NUMBER_OF_VARIABLES_QSS, 0.0);
  std::fill_n(current_substep.state_tq, NUMBER_OF_VARIABLES_QSS, 0.0);

  current_substep.t = 0;
  current_substep.extrapolation_method = qss_order;

  update_field();
  for (G4int i = 0; i < NUMBER_OF_VARIABLES_QSS; ++i)
  {
    dq_vector[i] = std::fmax(dqmin[INDEX_TYPE(i)], dqrel[INDEX_TYPE(i)] * std::fabs(current_substep.state_x[DERIVATIVE_0][i]));
    update_x_derivates_using_q(i);
    update_sync_time(i);
  }
}

// ----------------------------------------------------------------------------

void G4QSStepper::update_sync_time(G4int index)
{
  G4double &dq = dq_vector[index];
  G4double delta_sync_t = INFTY;
  // polynomial coefficients in increasing order of power, constant, linear, quadratic, etc
  G4double c, b, a, h;

  a = current_substep.state_x[DERIVATIVE_2][index]/2;
  b = current_substep.state_x[DERIVATIVE_1][index] - current_substep.state_q[DERIVATIVE_1][index] ;
  c = current_substep.state_x[DERIVATIVE_0][index] - current_substep.state_q[DERIVATIVE_0][index];

  // third order polynomial. It's a long algorithm but not a complex one
  if (qss_order == 3 && current_substep.state_x[DERIVATIVE_3][index] != 0.0)
  {
    // extra coefficient and h for the cubic polynomial and inclusion of the second order term from q
    h = current_substep.state_x[DERIVATIVE_3][index]/6;
    a -= current_substep.state_q[DERIVATIVE_2][index]/2;

    G4double q_cube = fCharge*fCharge*fCharge;
    
    // special case of | h * t3  | = dq
    if (a == 0 && b == 0 && c == 0)
    {
      delta_sync_t = cbrt(std::fabs(dq/h));
    }
    else
    {
      a /= h;
      b /= h;
      c /= h;

      G4double qLocal = (a * a - 3 * b) * (1.0 / 9.0);
      G4double r_base = (2*a*a*a - 9*a*b + 27*c)*(1.0/54.0);

      G4double sqrt_q = std::sqrt(qLocal);
      G4double sqrt_q_cube = std::sqrt(q_cube);
      G4double a_over_3 = a/3;
      for (G4double dQ : {dq,-dq})
      {
        G4double r = r_base + dQ/(2*h);
        // three real roots
        if (r*r < q_cube)
        {
          G4double theta = std::acos(r/sqrt_q_cube);
          G4double t1 = -2*sqrt_q*std::cos((1./3.)*theta) - a_over_3;
          G4double t2 = -2*sqrt_q*std::cos((1./3.)*(theta+2*CLHEP::pi)) - a_over_3;
          G4double t3 = -2*sqrt_q*std::cos((1./3.)*(theta-2*CLHEP::pi)) - a_over_3;
          for (G4double t : {t1,t2,t3})
          {
            if (t > 0) { delta_sync_t = std::fmin(delta_sync_t,t); }
          }
        }
        // one real root
        else
        {
          G4double A = -copysign(1,r) * cbrt(std::fabs(r) + std::sqrt(r*r - q_cube));
          G4double B = A == 0 ? 0 : qLocal/A;
          G4double t1 = A + B - a_over_3;
          if (t1 > 0) {delta_sync_t = std::fmin(delta_sync_t,t1);}
        }
      }
    }
  }

  // first order polynomial
  else if (qss_order == 1 || a == 0)
  {
    // dq = | b * t + c |
    if (b == 0) { delta_sync_t =  INFTY; }
    // (dq-c)/b > 0 <--> (b > 0 && dq > c) || (b < 0 && dq < c)
    // so we use dq if any of the cases holds and -dq if not
    else if ( (b > 0) == (dq > c) ) { delta_sync_t =  (dq-c)/b; }
    else { delta_sync_t = (-dq-c)/b; }
  }
  // second order polynomial
  else
  {
    if (b == 0)
    {
      // dq = | a_x * t2 + c |
      // identical to first order case but with sqrt
      if ((a > 0) == (dq > c)) {delta_sync_t =  std::sqrt((dq-c)/a);}
      else {delta_sync_t = std::sqrt((-dq-c)/a);}
    }
    else
    {
      // check both discriminants for both dq and - dq
      G4double a4 = 4*a;
      G4double a2 = 2*a;
      G4double discriminator_base = b*b - a4*c;
      G4double discriminator_difference = a4*dq;
      G4double discriminator_1 = discriminator_base + discriminator_difference;
      G4double discriminator_2 = discriminator_base - discriminator_difference;
      G4double fixed_solution_part = -b/a2;

      // simple trick to combine answers from all 4 solutions
      for(G4double discriminator : {discriminator_1, discriminator_2})
      {
        if (discriminator < 0) { continue; }
        G4double variable_solution_part = std::sqrt(discriminator)/std::fabs(a2);
        G4double t_local = fixed_solution_part - variable_solution_part;

        if (t_local <= 0 )
        {
          t_local = fixed_solution_part + variable_solution_part;
        }
        if (t_local > 0) { delta_sync_t = std::fmin(delta_sync_t,t_local); }
      }
    }
  }

  current_substep.sync_t[index] = current_substep.state_tx[index] + delta_sync_t;
}

// ----------------------------------------------------------------------------

void G4QSStepper::Stepper( const G4double y[],
	                         const G4double /*dydx*/ [],
	                         G4double h,
	                         G4double yout[],
	                         G4double /* yerr */ [] )
{
  using std::memcpy;
 
  initialize(y);

  const G4int QSS_MAX_SUBSTEPS = G4QSSMessenger::instance()->GetMaxSubsteps();

  G4double t = 0;

  fFinal_t = h/fVelocity;
  fFinal_t = std::fmin(fFinal_t,INFTY);

  while (t < fFinal_t && t < INFTY && substeps.current_substep_index < QSS_MAX_SUBSTEPS)
  {
    substeps.save_substep(&current_substep);

    // get minimum that makes some variable get too far from its quantized version
    G4int sync_index = get_next_sync_index();
    t = current_substep.sync_t[sync_index];
    t = std::fmin(t,fFinal_t);
    current_substep.t = t;

    // sync both and update their data
    // update x
    update_x(sync_index,t);

    // sync q
    current_substep.state_q[DERIVATIVE_0][sync_index] = current_substep.state_x[DERIVATIVE_0][sync_index];
    current_substep.state_q[DERIVATIVE_1][sync_index] = current_substep.state_x[DERIVATIVE_1][sync_index];
    current_substep.state_q[DERIVATIVE_2][sync_index] = current_substep.state_x[DERIVATIVE_2][sync_index];

    current_substep.state_tq[sync_index] = current_substep.state_tx[sync_index];

    dq_vector[sync_index] = std::fmax(dqmin[INDEX_TYPE(sync_index)], dqrel[INDEX_TYPE(sync_index)] * std::fabs(current_substep.state_x[DERIVATIVE_0][sync_index]));


    // Somehow this seems to be faster than the one below
    update_sync_time(sync_index);
    //  the trick belows work but seems to be slower
    //update_sync_time_one_coefficient(sync_index);


    // only update field if we actually changed position, not velocity
    // previous version called this every time which is unnecessary if field constant, and we bite the bullet if not
    if (sync_index < VELOCITY_IDX) { update_field(); }

    // we need to update the affected derivates of the other states
    G4double &tIndex = current_substep.state_tx[sync_index];


    // if we update position but magnetic field hasn't change then no other variables are affected!
    if(sync_index < VELOCITY_IDX && ! fField_changed) { continue; }

    // as qs are in different ts, we need to extrapolate the needed qs
    // we always need to extrapolate the velocity ones (because the lorentz equation)
    update_q(VX,tIndex);
    update_q(VY,tIndex);
    update_q(VZ,tIndex);


    // check which equations are altered by this update according to lorentz eq

    // b-field changed, need to update velocity states derivates
    if (sync_index < VELOCITY_IDX)
    {
      for (G4int i = VELOCITY_IDX; i < 6; ++i)
      {
        update_x(i,tIndex);
        update_x_velocity_derivates_using_q(i);
        update_sync_time(i);
      }
    }

    // velocity changed, need the other velocity states derivates and the corresponding position one
    else
    {
      G4int indexDep1 = (sync_index + 2)%VELOCITY_IDX  + VELOCITY_IDX;
      G4int indexDep2 = (sync_index + 1)%VELOCITY_IDX  + VELOCITY_IDX;
      G4int index_class = sync_index - VELOCITY_IDX;
      update_q(index_class,tIndex); // not updated before so we need to update it
      for (G4int i : {indexDep1, indexDep2, index_class})
      {
        update_x(i,tIndex);
        update_x_derivates_using_q(i);
        update_sync_time(i);
      }
    }
  }
  if(substeps.current_substep_index >= QSS_MAX_SUBSTEPS)
  {
    fFinal_t = current_substep.t;
  }

  for (G4int i = 0; i < NUMBER_OF_VARIABLES_QSS; ++i)
  {
    update_x(i, fFinal_t);
  }
  memcpy(yout, &current_substep.state_x[DERIVATIVE_0], sizeof(QSStateVector));

  velocity_to_momentum(yout);

  // fyout is used by interpolation driver, so we have to do this
  memcpy(fYout,yout,NUMBER_OF_VARIABLES_QSS*sizeof(G4double));
}

// ----------------------------------------------------------------------------

void G4QSStepper::Interpolate(G4double tau,G4double yOut[])
{
   G4double target_t = current_substep.t  * tau;
   G4int i = 0;
    G4double t = current_substep.t  * tau;;
  // linear search
  if (substeps.current_substep_index < 20)
  {
    while(i < substeps.current_substep_index && substeps._substeps[i+1].t  <= target_t )
    {
      i++;
    }
  }
  // binary search
  else 
  {
    G4int high_i = substeps.current_substep_index;
    G4int low_i = 0;
    G4int idx = high_i >> 1;
    while(low_i < high_i-1)
    {
      if(target_t < substeps._substeps[idx].t)
      {
        high_i = idx;
      }
      else
      {
        low_i = idx;
      }
      idx = (low_i+high_i) >> 1;
    }
    i = low_i;
  }

  extrapolate_all_states_to_t(&substeps._substeps[i], t, yOut);

  velocity_to_momentum(yOut);
}
