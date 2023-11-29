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
// G4QSS2
//
// G4QSS2 simulator

// Authors: Lucio Santi, Rodrigo Castro - 2018-2021
// --------------------------------------------------------------------
#ifndef _G4QSS2_H_
#define _G4QSS2_H_

#include "G4Types.hh"  //  For G4int, G4double
#include "G4qss_misc.hh"

#include <cmath>
#include <cassert>

#define  REPORT_CRITICAL_PROBLEM  1

#ifdef   REPORT_CRITICAL_PROBLEM
#include <cassert>
#include "G4Log.hh"
#endif

class G4QSS2
{
  public:

    G4QSS2(QSS_simulator sim) : simulator(sim) {}

    inline QSS_simulator getSimulator() const { return this->simulator; }

    inline G4int order() const { return 2; }

    inline void full_definition(G4double coeff)
    {
      G4double* const x = simulator->q;
      G4double* const dx = simulator->x;
      G4double* const alg = simulator->alg;

      dx[1] = x[9];
      dx[2] = 0;

      dx[4] = x[12];
      dx[5] = 0;

      dx[7] = x[15];
      dx[8] = 0;

      dx[10] = coeff * (alg[2] * x[12] - alg[1] * x[15]);
      dx[11] = 0;

      dx[13] = coeff * (alg[0] * x[15] - alg[2] * x[9]);
      dx[14] = 0;

      dx[16] = coeff * (alg[1] * x[9] - alg[0] * x[12]);
      dx[17] = 0;
    }

    inline void dependencies(G4int i, G4double coeff)
    {
      G4double* const x = simulator->q;
      G4double* const der = simulator->x;
      G4double* const alg = simulator->alg;

      switch (i)
      {
        case 0:
          der[10] = coeff * (alg[2] * x[12] - alg[1] * x[15]);
          der[11] = ((alg[2] * x[13] - x[16] * alg[1]) * coeff) / 2;

          der[13] = coeff * (alg[0] * x[15] - alg[2] * x[9]);
          der[14] = ((alg[0] * x[16] - alg[2] * x[10]) * coeff) / 2;

          der[16] = coeff * (alg[1] * x[9] - alg[0] * x[12]);
          der[17] = (-coeff * (alg[0] * x[13] - x[10] * alg[1])) / 2;
          return;
        case 1:
          der[10] = coeff * (alg[2] * x[12] - alg[1] * x[15]);
          der[11] = ((alg[2] * x[13] - x[16] * alg[1]) * coeff) / 2;

          der[13] = coeff * (alg[0] * x[15] - alg[2] * x[9]);
          der[14] = ((alg[0] * x[16] - alg[2] * x[10]) * coeff) / 2;

          der[16] = coeff * (alg[1] * x[9] - alg[0] * x[12]);
          der[17] = (-coeff * (alg[0] * x[13] - x[10] * alg[1])) / 2;
          return;
        case 2:
          der[10] = coeff * (alg[2] * x[12] - alg[1] * x[15]);
          der[11] = ((alg[2] * x[13] - x[16] * alg[1]) * coeff) / 2;

          der[13] = coeff * (alg[0] * x[15] - alg[2] * x[9]);
          der[14] = ((alg[0] * x[16] - alg[2] * x[10]) * coeff) / 2;

          der[16] = coeff * (alg[1] * x[9] - alg[0] * x[12]);
          der[17] = (-coeff * (alg[0] * x[13] - x[10] * alg[1])) / 2;
          return;
        case 3:
          der[1] = x[9];
          der[2] = (x[10]) / 2;

          der[13] = coeff * (alg[0] * x[15] - alg[2] * x[9]);
          der[14] = ((alg[0] * x[16] - alg[2] * x[10]) * coeff) / 2;

          der[16] = coeff * (alg[1] * x[9] - alg[0] * x[12]);
          der[17] = (-coeff * (alg[0] * x[13] - x[10] * alg[1])) / 2;
          return;
        case 4:
          der[4] = x[12];
          der[5] = (x[13]) / 2;

          der[10] = coeff * (alg[2] * x[12] - alg[1] * x[15]);
          der[11] = ((alg[2] * x[13] - x[16] * alg[1]) * coeff) / 2;

          der[16] = coeff * (alg[1] * x[9] - alg[0] * x[12]);
          der[17] = (-coeff * (alg[0] * x[13] - x[10] * alg[1])) / 2;
          return;
        case 5:
          der[7] = x[15];
          der[8] = (x[16]) / 2;

          der[10] = coeff * (alg[2] * x[12] - alg[1] * x[15]);
          der[11] = ((alg[2] * x[13] - x[16] * alg[1]) * coeff) / 2;

          der[13] = coeff * (alg[0] * x[15] - alg[2] * x[9]);
          der[14] = ((alg[0] * x[16] - alg[2] * x[10]) * coeff) / 2;
          return;
      }
    }

    inline void recompute_next_times(G4int* inf, G4double t)
    {
      G4int i;
      G4double* x = simulator->x;
      G4double* q = simulator->q;
      G4double* lqu = simulator->lqu;
      G4double* time = simulator->nextStateTime;

      for (i = 0; i < 3; i++)
      {
        const G4int var = inf[i];
        const G4int icf0 = 3 * var;
        const G4int icf1 = icf0 + 1;
        const G4int icf2 = icf1 + 1;

        time[var] = t;

        if (std::fabs(q[icf0] - x[icf0]) < lqu[var])
        {
          G4double mpr = -1, mpr2;
          G4double cf0 = q[icf0] + lqu[var] - x[icf0];
          G4double cf1 = q[icf1] - x[icf1];
          G4double cf2 = -x[icf2];
          G4double cf0Alt = q[icf0] - lqu[var] - x[icf0];

          if (unlikely(cf2 == 0 || (1000 * std::fabs(cf2)) < std::fabs(cf1)))
          {
            if (cf1 == 0) {
              mpr = Qss_misc::INF;
            } else
            {
              mpr = -cf0 / cf1;
              mpr2 = -cf0Alt / cf1;
              if (mpr < 0 || (mpr2 > 0 && mpr2 < mpr)) { mpr = mpr2; }
            }

            if (mpr < 0) { mpr = Qss_misc::INF; }
          }
          else
          {
            static G4ThreadLocal unsigned long long okCalls=0LL, badCalls= 0LL;
            constexpr G4double dangerZone = 1.0e+30;
            static G4ThreadLocal G4double bigCf1_pr = dangerZone,
                                          bigCf2_pr = dangerZone;
            static G4ThreadLocal G4double bigCf1 = 0.0, bigCf2 = 0.0;
            if( std::abs(cf1) > dangerZone || std::fabs(cf2) > dangerZone )
            {
              badCalls++;
              if( badCalls == 1 
                 || ( badCalls < 1000 && badCalls % 20 == 0 )
                 || (   1000 < badCalls && badCalls <   10000 && badCalls %  100 == 0 )
                 || (  10000 < badCalls && badCalls <  100000 && badCalls % 1000 == 0 )
                 || ( 100000 < badCalls &&                       badCalls % 10000 == 0 )
                 || ( std::fabs(cf1) > 1.5 * bigCf1_pr || std::fabs(cf2) > 1.5 * bigCf2_pr )
                )
              {
                std::cout << " cf1 = " << std::setw(15) << cf1 << " cf2= " << std::setw(15) << cf2
                          << "  badCall # " << badCalls << " of " << badCalls + okCalls
                          << "  fraction = " << double(badCalls) / double(badCalls+okCalls);

                if( std::fabs(cf1) > 1.5 * bigCf1_pr ) { bigCf1_pr = std::fabs(cf1); std::cout << " Bigger cf1 "; }
                if( std::fabs(cf2) > 1.5 * bigCf2_pr ) { bigCf2_pr = std::fabs(cf2); std::cout << " Bigger cf2 "; }
                std::cout << std::endl;
              }
              if( std::fabs(cf1) > 1.5 * bigCf1 ) { bigCf1 = std::fabs(cf1); }
              if( std::fabs(cf2) > 1.5 * bigCf2 ) { bigCf2 = std::fabs(cf2); }
            }
            else
            {
              okCalls++;
            }

#ifdef REPORT_CRITICAL_PROBLEM
            constexpr unsigned int exp_limit= 140;
            constexpr G4double limit= 1.0e+140; // std::pow(10,exp_limit));
            assert( std::fabs( std::pow(10, exp_limit) - limit ) < 1.0e-14*limit );
            G4bool bad_cf2fac= G4Log(std::fabs(cf2))
                             + G4Log(std::max( std::fabs(cf0), std::fabs(cf0Alt))) > 2*limit;
            if( std::fabs(cf1) > limit
               || G4Log(std::fabs(cf2))
                + G4Log(std::max( std::fabs(cf0), std::fabs(cf0Alt))) > 2*exp_limit )
            {
              G4ExceptionDescription ermsg;
              ermsg << "QSS2: Coefficients exceed tolerable values -- beyond " << limit << G4endl;
              if( std::fabs(cf1) > limit )
              {
                ermsg << " |cf1| = " << cf1 << " is > " << limit << " (limit)";
              }
              if( bad_cf2fac)
              {
                ermsg << " bad cf2-factor:  cf2 = " << cf2
                      << " product is > " << 2*limit << " (limit)";
              }
              G4Exception("QSS2::recompute_next_times",
                          "Field/Qss2-", FatalException, ermsg ); 
            }
#endif
            G4double cf1_2 = cf1 * cf1;
            G4double cf2_4 = 4 * cf2;
            G4double disc1 = cf1_2 - cf2_4 * cf0;
            G4double disc2 = cf1_2 - cf2_4 * cf0Alt;
            G4double cf2_d2 = 2 * cf2;

            if (unlikely(disc1 < 0 && disc2 < 0))  // no real roots
            {
              mpr = Qss_misc::INF;
            }
            else if (disc2 < 0)
            {
              G4double sd, r1;
              sd = std::sqrt(disc1);
              r1 = (-cf1 + sd) / cf2_d2;
              if (r1 > 0) {
                mpr = r1;
              } else {
                mpr = Qss_misc::INF;
              }
              r1 = (-cf1 - sd) / cf2_d2;
              if ((r1 > 0) && (r1 < mpr)) { mpr = r1; }
            }
            else if (disc1 < 0)
            {
              G4double sd, r1;
              sd = std::sqrt(disc2);
              r1 = (-cf1 + sd) / cf2_d2;
              if (r1 > 0) {
                mpr = r1;
              } else {
                mpr = Qss_misc::INF;
              }
              r1 = (-cf1 - sd) / cf2_d2;
              if ((r1 > 0) && (r1 < mpr)) { mpr = r1; }
            }
            else
            {
              G4double sd1, r1, sd2, r2;
              sd1 = std::sqrt(disc1);
              sd2 = std::sqrt(disc2);
              r1 = (-cf1 + sd1) / cf2_d2;
              r2 = (-cf1 + sd2) / cf2_d2;
              if (r1 > 0) { mpr = r1; }
              else { mpr = Qss_misc::INF; }
              r1 = (-cf1 - sd1) / cf2_d2;
              if ((r1 > 0) && (r1 < mpr)) { mpr = r1; }
              if (r2 > 0 && r2 < mpr) { mpr = r2; }
              r2 = (-cf1 - sd2) / cf2_d2;
              if ((r2 > 0) && (r2 < mpr)) { mpr = r2; }
            }
          }
          time[var] += mpr;
        }
      }
    }

    inline void recompute_all_state_times(G4double t)
    {
      G4double mpr;
      G4double* const x = simulator->x;
      G4double* const lqu = simulator->lqu;
      G4double* const time = simulator->nextStateTime;

      for (G4int var = 0, icf0 = 0; var < 6; var++, icf0 += 3)
      {
        const G4int icf1 = icf0 + 1;

        if (x[icf1] == 0)
        {
          time[var] = Qss_misc::INF;
        }
        else
        {
          mpr = lqu[var] / x[icf1];
          if (mpr < 0) { mpr *= -1; }
          time[var] = t + mpr;
        }
      }
    }

    inline void next_time(G4int var, G4double t)
    {
      const G4int cf2 = var * 3 + 2;
      G4double* const x = simulator->x;
      G4double* const lqu = simulator->lqu;
      G4double* const time = simulator->nextStateTime;

      if (x[cf2] != 0.0) {
        time[var] = t + std::sqrt(lqu[var] / std::fabs(x[cf2]));
      } else {
        time[var] = Qss_misc::INF;
      }
    }

    inline void update_quantized_state(G4int i)
    {
      const G4int cf0 = i * 3, cf1 = cf0 + 1;
      G4double* const q = simulator->q;
      G4double* const x = simulator->x;

      q[cf0] = x[cf0];
      q[cf1] = x[cf1];
    }

    inline void reset_state(G4int i, G4double value)
    {
      G4double* const x = simulator->x;
      G4double* const q = simulator->q;
      G4double* const tq = simulator->tq;
      G4double* const tx = simulator->tx;
      const G4int idx = 3 * i;

      x[idx] = value;

      simulator->lqu[i] = simulator->dQRel[i] * std::fabs(value);
      if (simulator->lqu[i] < simulator->dQMin[i])
      {
        simulator->lqu[i] = simulator->dQMin[i];
      }

      q[idx] = value;
      q[idx + 1] = tq[i] = tx[i] = 0;
    }

    inline G4double evaluate_x_poly(G4int i, G4double dt, G4double* p)
    {
      return (p[i + 2] * dt + p[i + 1]) * dt + p[i];
    }

    inline void advance_time_q(G4int i, G4double dt)  // __attribute__((hot))
    {
      G4double* const p = simulator->q;
      p[i] = p[i] + dt * p[i + 1];
    }

    inline void advance_time_x(G4int i, G4double dt)  // __attribute__((hot))
    {
      G4double* const p = simulator->x;
      const G4int i0 = i, i1 = i0 + 1, i2 = i1 + 1;
      p[i0] = (p[i2] * dt + p[i1]) * dt + p[i0];
      p[i1] = p[i1] + 2 * dt * p[i2];
    }

  private:

    QSS_simulator simulator;
};

#endif
