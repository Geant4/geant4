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
// G4BulirschStoer class implementation
// Based on bulirsch_stoer.hpp from boost
//
// Author: Dmitry Sorokin, Google Summer of Code 2016
// --------------------------------------------------------------------

#include "G4BulirschStoer.hh"

#include "G4FieldUtils.hh"

namespace
{
  constexpr G4double STEPFAC1 = 0.65;
  constexpr G4double STEPFAC2 = 0.94;
  constexpr G4double STEPFAC3 = 0.02;
  constexpr G4double STEPFAC4 = 4.0;
  constexpr G4double KFAC1 = 0.8;
  constexpr G4double KFAC2 = 0.9;
  constexpr G4double inv_STEPFAC1 = 1.0 / STEPFAC1;
  constexpr G4double inv_STEPFAC4 = 1.0 / STEPFAC4;
} // namespace

G4BulirschStoer::G4BulirschStoer(G4EquationOfMotion* equation,
                                 G4int nvar, G4double eps_rel, G4double max_dt)
  : fnvar(nvar), m_eps_rel(eps_rel), m_midpoint(equation,nvar),
    m_last_step_rejected(false), m_first(true), m_dt_last(0.0), m_max_dt(max_dt)
{
  /* initialize sequence of stage numbers and work */

  for(G4int i = 0; i < m_k_max + 1; ++i)
  {
    m_interval_sequence[i] = 2 * (i + 1);
    if (i == 0)
    {
      m_cost[i] = m_interval_sequence[i];
    }
    else
    {
      m_cost[i] = m_cost[i-1] + m_interval_sequence[i];
    }
    for(G4int k = 0; k < i; ++k)
    {
      const G4double r = static_cast<G4double>(m_interval_sequence[i])
                       / static_cast<G4double>(m_interval_sequence[k]);
      m_coeff[i][k] = 1.0 / (r * r - 1.0); // coefficients for extrapolation
    }

    // crude estimate of optimal order
    m_current_k_opt = 4;

    // no calculation because log10 might not exist for value_type!

    //const G4double logfact = -log10(std::max(eps_rel, 1.0e-12)) * 0.6 + 0.5;
    //m_current_k_opt = std::max(1.,
    //                  std::min(static_cast<G4double>(m_k_max-1), logfact));
  }
}

G4BulirschStoer::step_result
G4BulirschStoer::try_step( const G4double in[], const G4double dxdt[],
                           G4double& t, G4double out[], G4double& dt)
{
  if(m_max_dt < dt)
  {
    // given step size is bigger then max_dt set limit and return fail
    //
    dt = m_max_dt;
    return step_result::fail;
  }

  if (dt != m_dt_last)
  {
    reset(); // step size changed from outside -> reset
  }

  G4bool reject = true;

  G4double new_h = dt;

  /* m_current_k_opt is the estimated current optimal stage number */

  for(G4int k = 0; k <= m_current_k_opt+1; ++k)
  {
    // the stage counts are stored in m_interval_sequence
    //
    m_midpoint.SetSteps(m_interval_sequence[k]);
    if(k == 0)
    {
      m_midpoint.DoStep(in, dxdt, out, dt);
      /* the first step, nothing more to do */
    }
    else
    {
      m_midpoint.DoStep(in, dxdt, m_table[k-1], dt);
      extrapolate(k, out);
      // get error estimate
      for (G4int i = 0; i < fnvar; ++i)
      {
        m_err[i] = out[i] - m_table[0][i];
      }
      const G4double error =
            field_utils::relativeError(out, m_err, dt, m_eps_rel);
      h_opt[k] = calc_h_opt(dt, error, k);
      work[k] = static_cast<G4double>(m_cost[k]) / h_opt[k];

      if( (k == m_current_k_opt-1) || m_first)  // convergence before k_opt ?
      {
        if(error < 1.0)
        {
          // convergence
          reject = false;
          if( (work[k] < KFAC2 * work[k-1]) || (m_current_k_opt <= 2) )
          {
            // leave order as is (except we were in first round)
            m_current_k_opt = std::min(m_k_max - 1 , std::max(2 , k + 1));
            new_h = h_opt[k];
            new_h *= static_cast<G4double>(m_cost[k + 1])
                   / static_cast<G4double>(m_cost[k]);
          }
          else
          {
            m_current_k_opt = std::min(m_k_max - 1, std::max(2, k));
            new_h = h_opt[k];
          }
          break;
        }
        else if(should_reject(error , k) && !m_first)
        {
          reject = true;
          new_h = h_opt[k];
          break;
        }
      }
      if(k == m_current_k_opt)  // convergence at k_opt ?
      {
        if(error < 1.0)
        {
          // convergence
          reject = false;
          if(work[k-1] < KFAC2 * work[k])
          {
            m_current_k_opt = std::max( 2 , m_current_k_opt-1 );
            new_h = h_opt[m_current_k_opt];
          }
          else if( (work[k] < KFAC2 * work[k-1]) && !m_last_step_rejected )
          {
            m_current_k_opt = std::min(m_k_max - 1, m_current_k_opt + 1);
            new_h = h_opt[k];
            new_h *= static_cast<G4double>(m_cost[m_current_k_opt])
                   / static_cast<G4double>(m_cost[k]);
          }
          else
          {
            new_h = h_opt[m_current_k_opt];
          }
          break;
        }
        else if(should_reject(error, k))
        {
          reject = true;
          new_h = h_opt[m_current_k_opt];
          break;
        }
      }
      if(k == m_current_k_opt + 1)  // convergence at k_opt+1 ?
      {
        if(error < 1.0)  // convergence
        {
          reject = false;
          if(work[k-2] < KFAC2 * work[k-1])
          {
            m_current_k_opt = std::max(2, m_current_k_opt - 1);
          }
          if((work[k] < KFAC2 * work[m_current_k_opt]) && !m_last_step_rejected)
          {
            m_current_k_opt = std::min(m_k_max - 1 , k);
          }
          new_h = h_opt[m_current_k_opt];
        }
        else
        {
          reject = true;
          new_h = h_opt[m_current_k_opt];
        }
        break;
      }
    }
  }

  if(!reject)
  {
    t += dt;
  }

  if(!m_last_step_rejected || new_h < dt)
  {
    // limit step size
    new_h = std::min(m_max_dt, new_h);
    m_dt_last = new_h;
    dt = new_h;
  }

  m_last_step_rejected = reject;
  m_first = false;

  return reject ? step_result::fail : step_result::success;
}

void G4BulirschStoer::reset()
{
  m_first = true;
  m_last_step_rejected = false;
}

void G4BulirschStoer::extrapolate(std::size_t k , G4double xest[])
{
  /* polynomial extrapolation, see http://www.nr.com/webnotes/nr3web21.pdf
   * uses the obtained intermediate results to extrapolate to dt->0 */

  for(std::size_t j = k - 1 ; j > 0; --j)
  {
    for (G4int i = 0; i < fnvar; ++i)
    {
      m_table[j-1][i] = m_table[j][i] * (1. + m_coeff[k][j])
                      - m_table[j-1][i] * m_coeff[k][j];
    }
  }
  for (G4int i = 0; i < fnvar; ++i)
  {
    xest[i] = m_table[0][i] * (1. + m_coeff[k][0]) - xest[i] * m_coeff[k][0];
  }
}

G4double
G4BulirschStoer::calc_h_opt(G4double h , G4double error , std::size_t k) const
{
  /* calculates the optimal step size for a given error and stage number */

  const G4double expo =  1.0 / (2 * k + 1);
  const G4double facmin = std::pow(STEPFAC3, expo);
  G4double fac;

  G4double facminInv= 1.0 / facmin;
  if (error == 0.0)
  {
    fac = facminInv;
  }
  else
  {
    fac = STEPFAC2 * std::pow(error * inv_STEPFAC1 , -expo);
    fac = std::max(facmin * inv_STEPFAC4, std::min( facminInv, fac));
  }

  return h * fac;
}

//why is not used!!??
G4bool G4BulirschStoer::set_k_opt(std::size_t k, G4double& dt)
{
  /* calculates the optimal stage number */

  if(k == 1)
  {
    m_current_k_opt = 2;
    return true;
  }
  if( (work[k-1] < KFAC1 * work[k]) || (k == m_k_max) )   // order decrease
  {
    m_current_k_opt = (G4int)k - 1;
    dt = h_opt[ m_current_k_opt ];
    return true;
  }
  else if( (work[k] < KFAC2 * work[k-1])
          || m_last_step_rejected || (k == m_k_max-1) )
  {  // same order - also do this if last step got rejected
    m_current_k_opt = (G4int)k;
    dt = h_opt[m_current_k_opt];
    return true;
  }
  else {   // order increase - only if last step was not rejected
    m_current_k_opt = (G4int)k + 1;
    dt = h_opt[m_current_k_opt - 1] * m_cost[m_current_k_opt]
       / m_cost[m_current_k_opt - 1];
    return true;
  }
}

G4bool G4BulirschStoer::in_convergence_window(G4int k) const
{
  if( (k == m_current_k_opt - 1) && !m_last_step_rejected )
  {
    return true; // decrease stepsize only if last step was not rejected
  }
  return (k == m_current_k_opt) || (k == m_current_k_opt + 1);
}


G4bool G4BulirschStoer::should_reject(G4double error, G4int k) const
{
  if(k == m_current_k_opt - 1)
  {
    const G4double d = G4double(m_interval_sequence[m_current_k_opt]
                              * m_interval_sequence[m_current_k_opt+1]);
    const G4double e = G4double(m_interval_sequence[0]);
    const G4double e2 = e*e; 
    // step will fail, criterion 17.3.17 in NR
    return error * e2 * e2 > d * d;  //  was return error > dOld * dOld; (where dOld= d/e; )
  }
  else if(k == m_current_k_opt)
  {
    const G4double d = G4double(m_interval_sequence[m_current_k_opt]);
    const G4double e = G4double(m_interval_sequence[0]);
    return error * e * e > d * d; //  was return error > dOld * dOld; (where dOld= d/e; )
  }
  else
  {
    return error > 1.0;
  }
}
