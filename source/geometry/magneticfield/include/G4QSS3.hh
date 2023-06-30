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
// G4QSS3
//
// G4QSS3 simulator

// Authors: Lucio Santi, Rodrigo Castro - 2018-2021
// --------------------------------------------------------------------
#ifndef _G4QSS3_H_
#define _G4QSS3_H_

#include "G4Types.hh"
#include "G4qss_misc.hh"

#include <cmath>

class G4QSS3
{
  public:

    G4QSS3(QSS_simulator);

    inline QSS_simulator getSimulator() const { return this->simulator; }

    inline G4int order() const { return 3; }

    inline void full_definition(G4double coeff)
    {
      G4double* const x = simulator->q;
      G4double* const dx = simulator->x;
      G4double* const alg = simulator->alg;

      dx[1] = x[12];
      dx[2] = 0;
      dx[3] = 0;

      dx[5] = x[16];
      dx[6] = 0;
      dx[7] = 0;

      dx[9] = x[20];
      dx[10] = 0;
      dx[11] = 0;

      dx[13] = coeff * (alg[2] * x[16] - alg[1] * x[20]);
      dx[14] = 0;
      dx[15] = 0;

      dx[17] = coeff * (alg[0] * x[20] - alg[2] * x[12]);
      dx[18] = 0;
      dx[19] = 0;

      dx[21] = coeff * (alg[1] * x[12] - alg[0] * x[16]);
      dx[22] = 0;
      dx[23] = 0;
    }

    inline void dependencies(G4int i, G4double coeff)
    {
      G4double* const x = simulator->q;
      G4double* const der = simulator->x;
      G4double* const alg = simulator->alg;

      switch (i)
      {
        case 0:
          der[13] = coeff * (alg[2] * x[16] - alg[1] * x[20]);
          der[14] = ((alg[2] * x[17] - x[21] * alg[1]) * coeff) / 2;
          der[15] = (coeff * (alg[2] * x[18] - x[22] * alg[1])) / 3;

          der[17] = coeff * (alg[0] * x[20] - alg[2] * x[12]);
          der[18] = ((alg[0] * x[21] - alg[2] * x[13]) * coeff) / 2;
          der[19] = (coeff * (alg[0] * x[22] - alg[2] * x[14])) / 3;

          der[21] = coeff * (alg[1] * x[12] - alg[0] * x[16]);
          der[22] = (coeff * (x[13] * alg[1] - alg[0] * x[17])) / 2;
          der[23] = (coeff * (alg[1] * x[14] - x[18] * alg[0])) / 3;
          return;
        case 1:
          der[13] = coeff * (alg[2] * x[16] - alg[1] * x[20]);
          der[14] = ((alg[2] * x[17] - x[21] * alg[1]) * coeff) / 2;
          der[15] = (coeff * (alg[2] * x[18] - x[22] * alg[1])) / 3;

          der[17] = coeff * (alg[0] * x[20] - alg[2] * x[12]);
          der[18] = ((alg[0] * x[21] - alg[2] * x[13]) * coeff) / 2;
          der[19] = (coeff * (alg[0] * x[22] - alg[2] * x[14])) / 3;

          der[21] = coeff * (alg[1] * x[12] - alg[0] * x[16]);
          der[22] = (coeff * (x[13] * alg[1] - alg[0] * x[17])) / 2;
          der[23] = (coeff * (alg[1] * x[14] - x[18] * alg[0])) / 3;
          return;
        case 2:
          der[13] = coeff * (alg[2] * x[16] - alg[1] * x[20]);
          der[14] = ((alg[2] * x[17] - x[21] * alg[1]) * coeff) / 2;
          der[15] = (coeff * (alg[2] * x[18] - x[22] * alg[1])) / 3;

          der[17] = coeff * (alg[0] * x[20] - alg[2] * x[12]);
          der[18] = ((alg[0] * x[21] - alg[2] * x[13]) * coeff) / 2;
          der[19] = (coeff * (alg[0] * x[22] - alg[2] * x[14])) / 3;

          der[21] = coeff * (alg[1] * x[12] - alg[0] * x[16]);
          der[22] = (coeff * (x[13] * alg[1] - alg[0] * x[17])) / 2;
          der[23] = (coeff * (alg[1] * x[14] - x[18] * alg[0])) / 3;
          return;
        case 3:
          der[1] = x[12];
          der[2] = x[13] / 2;
          der[3] = x[14] / 3;

          der[17] = coeff * (alg[0] * x[20] - alg[2] * x[12]);
          der[18] = ((alg[0] * x[21] - alg[2] * x[13]) * coeff) / 2;
          der[19] = (coeff * (alg[0] * x[22] - alg[2] * x[14])) / 3;

          der[21] = coeff * (alg[1] * x[12] - alg[0] * x[16]);
          der[22] = (coeff * (x[13] * alg[1] - alg[0] * x[17])) / 2;
          der[23] = (coeff * (alg[1] * x[14] - x[18] * alg[0])) / 3;
          return;
        case 4:
          der[5] = x[16];
          der[6] = x[17] / 2;
          der[7] = x[18] / 3;

          der[13] = coeff * (alg[2] * x[16] - alg[1] * x[20]);
          der[14] = ((alg[2] * x[17] - x[21] * alg[1]) * coeff) / 2;
          der[15] = (coeff * (alg[2] * x[18] - x[22] * alg[1])) / 3;

          der[21] = coeff * (alg[1] * x[12] - alg[0] * x[16]);
          der[22] = (coeff * (x[13] * alg[1] - alg[0] * x[17])) / 2;
          der[23] = (coeff * (alg[1] * x[14] - x[18] * alg[0])) / 3;
          return;
        case 5:
          der[9] = x[20];
          der[10] = x[21] / 2;
          der[11] = x[22] / 3;

          der[13] = coeff * (alg[2] * x[16] - alg[1] * x[20]);
          der[14] = ((alg[2] * x[17] - x[21] * alg[1]) * coeff) / 2;
          der[15] = (coeff * (alg[2] * x[18] - x[22] * alg[1])) / 3;

          der[17] = coeff * (alg[0] * x[20] - alg[2] * x[12]);
          der[18] = ((alg[0] * x[21] - alg[2] * x[13]) * coeff) / 2;
          der[19] = (coeff * (alg[0] * x[22] - alg[2] * x[14])) / 3;
          return;
      }
    }

    void recompute_next_times(G4int* inf, G4double t);  // __attribute__((hot));

    inline void recompute_all_state_times(G4double t)
    {
      G4double mpr;
      G4double* const x = simulator->x;
      G4double* const lqu = simulator->lqu;
      G4double* const time = simulator->nextStateTime;

      for (G4int var = 0, icf0 = 0; var < 6; var++, icf0 += 4)
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

    inline void next_time(G4int i, G4double t)
    {
      const G4int cf3 = 4 * i + 3;
      G4double* const x = simulator->x;
      G4double* const lqu = simulator->lqu;
      G4double* const time = simulator->nextStateTime;

      if (likely(x[cf3])) {
        time[i] = t + std::cbrt(lqu[i] / std::fabs(x[cf3]));
      } else {
        time[i] = Qss_misc::INF;
      }
    }

    inline void update_quantized_state(G4int i)
    {
      const G4int cf0 = i * 4, cf1 = cf0 + 1, cf2 = cf1 + 1;
      G4double* const q = simulator->q;
      G4double* const x = simulator->x;

      q[cf0] = x[cf0];
      q[cf1] = x[cf1];
      q[cf2] = x[cf2];
    }

    inline void reset_state(G4int i, G4double value)
    {
      G4double* const x = simulator->x;
      G4double* const q = simulator->q;
      G4double* const tq = simulator->tq;
      G4double* const tx = simulator->tx;
      const G4int idx = 4 * i;

      x[idx] = value;

      simulator->lqu[i] = simulator->dQRel[i] * std::fabs(value);
      if (simulator->lqu[i] < simulator->dQMin[i])
      {
        simulator->lqu[i] = simulator->dQMin[i];
      }
      q[idx] = value;
      q[idx + 1] = q[idx + 2] = tq[i] = tx[i] = 0;
    }

    inline G4double evaluate_x_poly(G4int i, G4double dt, G4double* p)
    {
      return ((p[i + 3] * dt + p[i + 2]) * dt + p[i + 1]) * dt + p[i];
    }

    inline void advance_time_q(G4int i, G4double dt)  //  __attribute__((hot))
    {
      G4double* const p = simulator->q;
      const G4int i0 = i, i1 = i0 + 1, i2 = i1 + 1;
      p[i0] = (p[i2] * dt + p[i1]) * dt + p[i0];
      p[i1] = p[i1] + 2 * dt * p[i2];
    }

    inline void advance_time_x(G4int i, G4double dt)  // __attribute__((hot))
    {
      G4double* const p = simulator->x;
      const G4int i0 = i, i1 = i0 + 1, i2 = i1 + 1, i3 = i2 + 1;
      p[i0] = ((p[i3] * dt + p[i2]) * dt + p[i1]) * dt + p[i0];
      p[i1] = (3 * p[i3] * dt + 2 * p[i2]) * dt + p[i1];
      p[i2] = p[i2] + 3 * dt * p[i3];
    }

    G4double min_pos_root(G4double* coeff, G4int order);

    inline G4double min_pos_root_2(G4double* coeff)
    {
      G4double mpr = Qss_misc::INF;

      if (coeff[2] == 0 || (1000 * std::fabs(coeff[2])) < std::fabs(coeff[1]))
      {
        if (coeff[1] == 0) {
          mpr = Qss_misc::INF;
        } else {
          mpr = -coeff[0] / coeff[1];
        }

        if (mpr < 0) { mpr = Qss_misc::INF; }
      }
      else
      {
        G4double disc;
        disc = coeff[1] * coeff[1] - 4 * coeff[2] * coeff[0];
        if (disc < 0)   // no real roots
        {
          mpr = Qss_misc::INF;
        }
        else
        {
          G4double sd, r1;
          G4double cf2_d2 = 2 * coeff[2];

          sd = std::sqrt(disc);
          r1 = (-coeff[1] + sd) / cf2_d2;
          if (r1 > 0) {
            mpr = r1;
          } else {
            mpr = Qss_misc::INF; 
          }
          r1 = (-coeff[1] - sd) / cf2_d2;
          if ((r1 > 0) && (r1 < mpr)) { mpr = r1; }
        }
      }

      return mpr;
    }  // __attribute__((hot))

    inline G4double min_pos_root_3(G4double* coeff)
    {
      G4double mpr = Qss_misc::INF;
      static const G4double sqrt3 = std::sqrt(3);

      if ((coeff[3] == 0) || (1000 * std::fabs(coeff[3]) < std::fabs(coeff[2])))
      {
        mpr = min_pos_root_2(coeff);
      }
      else if (coeff[0] == 0)
      {
        if (coeff[1] == 0)
        {
          mpr = -coeff[2] / coeff[3];
        }
        else
        {
          coeff[0] = coeff[1];
          coeff[1] = coeff[2];
          coeff[2] = coeff[3];
          mpr = min_pos_root_2(coeff);
        }
      }
      else
      {
        G4double q, r, disc, q3;
        G4double val = coeff[2] / 3 / coeff[3];
        G4double cf32 = coeff[3] * coeff[3];
        G4double cf22 = coeff[2] * coeff[2];
        G4double denq = 9 * cf32;
        G4double denr = 6 * coeff[3] * denq;
        G4double rcomm = 9 * coeff[3] * coeff[2] * coeff[1] - 2 * cf22 * coeff[2];

        q = (3 * coeff[3] * coeff[1] - cf22) / denq;
        q3 = q * q * q;

        r = (rcomm - 27 * cf32 * coeff[0]) / denr;
        disc = q3 + r * r;
        mpr = Qss_misc::INF;

        if (disc >= 0)
        {
          G4double sd, sx, t, r1, rsd;
          sd = std::sqrt(disc);
          rsd = r + sd;
          if (rsd > 0) {
            sx = std::cbrt(rsd);
          } else {
            sx = -std::cbrt(std::fabs(rsd));
          }

          rsd = r - sd;
          if (rsd > 0) {
            t = std::cbrt(rsd);
          } else {
            t = -std::cbrt(std::fabs(rsd));
          }

          r1 = sx + t - val;

          if (r1 > 0) { mpr = r1; }
        }
        else
        {
          // three real roots
          G4double rho, th, rho13, costh3, sinth3, spt, smti32, r1;
          rho = std::sqrt(-q3);
          th = std::acos(r / rho);
          rho13 = std::cbrt(rho);
          costh3 = std::cos(th / 3);
          sinth3 = std::sin(th / 3);
          spt = rho13 * 2 * costh3;
          smti32 = -rho13 * sinth3 * sqrt3;
          r1 = spt - val;
          if (r1 > 0) { mpr = r1; }
          r1 = -spt / 2 - val + smti32;
          if ((r1 > 0) && (r1 < mpr)) { mpr = r1; }
          r1 = r1 - 2 * smti32;
          if ((r1 > 0) && (r1 < mpr)) { mpr = r1; }
        }
      }

      return mpr;
    }  // __attribute__((hot))

    inline G4double min_pos_root_2_alt(G4double* coeff, G4double cf0Alt)
    {
      G4double mpr = Qss_misc::INF;
      G4double mpr2;

      if (coeff[2] == 0 || (1000 * std::fabs(coeff[2])) < std::fabs(coeff[1]))
      {
        if (coeff[1] == 0)
        {
          mpr = Qss_misc::INF;
        }
        else
        {
          mpr = -coeff[0] / coeff[1];
          mpr2 = -cf0Alt / coeff[1];
          if (mpr < 0 || (mpr2 > 0 && mpr2 < mpr)) { mpr = mpr2; }
        }

        if (mpr < 0) { mpr = Qss_misc::INF; }
      }
      else
      {
        G4double cf1_2 = coeff[1] * coeff[1];
        G4double cf2_4 = 4 * coeff[2];
        G4double disc1 = cf1_2 - cf2_4 * coeff[0];
        G4double disc2 = cf1_2 - cf2_4 * cf0Alt;
        G4double cf2_d2 = 2 * coeff[2];

        if (unlikely(disc1 < 0 && disc2 < 0))
        {
          mpr = Qss_misc::INF;
        }
        else if (disc2 < 0)
        {
          G4double sd, r1;
          sd = std::sqrt(disc1);
          r1 = (-coeff[1] + sd) / cf2_d2;
          if (r1 > 0) {
            mpr = r1;
          } else {
            mpr = Qss_misc::INF;
          }
          r1 = (-coeff[1] - sd) / cf2_d2;
          if ((r1 > 0) && (r1 < mpr)) { mpr = r1; }
        }
        else if (disc1 < 0)
        {
          G4double sd, r1;
          sd = std::sqrt(disc2);
          r1 = (-coeff[1] + sd) / cf2_d2;
          if (r1 > 0) {
            mpr = r1;
          } else {
            mpr = Qss_misc::INF;
          }
          r1 = (-coeff[1] - sd) / cf2_d2;
          if ((r1 > 0) && (r1 < mpr)) { mpr = r1; }
        }
        else
        {
          G4double sd1, r1, sd2, r2;
          sd1 = std::sqrt(disc1);
          sd2 = std::sqrt(disc2);
          r1 = (-coeff[1] + sd1) / cf2_d2;
          r2 = (-coeff[1] + sd2) / cf2_d2;

          if (r1 > 0) {
            mpr = r1;
          } else {
            mpr = Qss_misc::INF;
          }
          r1 = (-coeff[1] - sd1) / cf2_d2;
          if ((r1 > 0) && (r1 < mpr)) { mpr = r1; }

          if (r2 > 0 && r2 < mpr) { mpr = r2; }
          r2 = (-coeff[1] - sd2) / cf2_d2;
          if ((r2 > 0) && (r2 < mpr)) { mpr = r2; }
        }
      }

      return mpr;
    }  // __attribute__((hot))

    inline G4double min_pos_root_3_alt(G4double* coeff, G4double cf0Alt)
    {
      G4double mpr = Qss_misc::INF;
      static const G4double sqrt3 = std::sqrt(3);

      if ((coeff[3] == 0) || (1000 * std::fabs(coeff[3]) < std::fabs(coeff[2])))
      {
        mpr = min_pos_root_2_alt(coeff, cf0Alt);
      }
      else if (coeff[0] == 0)
      {
        G4double mpr2;
        coeff[0] = cf0Alt;
        mpr = min_pos_root_3(coeff);

        if (coeff[1] == 0)
        {
          mpr2 = -coeff[2] / coeff[3];
        }
        else
        {
          coeff[0] = coeff[1];
          coeff[1] = coeff[2];
          coeff[2] = coeff[3];
          mpr2 = min_pos_root_2(coeff);
        }

        if (mpr2 > 0 && mpr2 < mpr) { mpr = mpr2; }
      }
      else if (cf0Alt == 0)
      {
        G4double mpr2;
        mpr = min_pos_root_3(coeff);

        if (coeff[1] == 0)
        {
          mpr2 = -coeff[2] / coeff[3];
        }
        else
        {
          coeff[0] = coeff[1];
          coeff[1] = coeff[2];
          coeff[2] = coeff[3];
          mpr2 = min_pos_root_2(coeff);
        }

        if (mpr2 > 0 && mpr2 < mpr) { mpr = mpr2; }
      }
      else
      {
        G4double q, r, rAlt, disc, discAlt, q3;
        G4double val = coeff[2] / 3 / coeff[3];
        G4double cf32 = coeff[3] * coeff[3];
        G4double cf22 = coeff[2] * coeff[2];
        G4double denq = 9 * cf32;
        G4double denr = 6 * coeff[3] * denq;
        G4double rcomm = 9 * coeff[3] * coeff[2] * coeff[1] - 2 * cf22 * coeff[2];

        q = (3 * coeff[3] * coeff[1] - cf22) / denq;
        q3 = q * q * q;

        r = (rcomm - 27 * cf32 * coeff[0]) / denr;
        rAlt = (rcomm - 27 * cf32 * cf0Alt) / denr;

        disc = q3 + r * r;
        discAlt = q3 + rAlt * rAlt;
        mpr = Qss_misc::INF;

        if (disc >= 0)
        {
          G4double sd, sx, t, r1, rsd;
          sd = std::sqrt(disc);
          rsd = r + sd;
          if (rsd > 0) {
            sx = std::cbrt(rsd);
          } else {
            sx = -std::cbrt(std::fabs(rsd));
          }

          rsd = r - sd;
          if (rsd > 0) {
            t = std::cbrt(rsd);
          } else {
            t = -std::cbrt(std::fabs(rsd));
          }

          r1 = sx + t - val;

          if (r1 > 0) { mpr = r1; }

          if (discAlt >= 0)
          {
            G4double sdAlt, sAlt, tAlt, r1Alt, rsdAlt;
            sdAlt = std::sqrt(discAlt);
            rsdAlt = rAlt + sdAlt;
            if (rsdAlt > 0) {
              sAlt = std::cbrt(rsdAlt);
            } else {
              sAlt = -std::cbrt(std::fabs(rsdAlt));
            }

            rsdAlt = rAlt - sdAlt;
            if (rsdAlt > 0) {
              tAlt = std::cbrt(rsdAlt);
            } else {
              tAlt = -std::cbrt(std::fabs(rsdAlt));
            }

            r1Alt = sAlt + tAlt - val;

            if (r1Alt > 0 && r1Alt < mpr) { mpr = r1Alt; }
          }
          else
          {
            G4double rho, th, rho13, costh3, sinth3, spt, smti32, r1Alt;

            rho = std::sqrt(-q3);
            th = std::acos(rAlt / rho);
            rho13 = std::cbrt(rho);
            costh3 = std::cos(th / 3);
            sinth3 = std::sin(th / 3);
            spt = rho13 * 2 * costh3;
            smti32 = -rho13 * sinth3 * sqrt3;
            r1Alt = spt - val;
            if (r1Alt > 0 && r1Alt < mpr) { mpr = r1Alt; }
            r1Alt = -spt / 2 - val + smti32;
            if (r1Alt > 0 && r1Alt < mpr) { mpr = r1Alt; }
            r1Alt = r1Alt - 2 * smti32;
            if (r1Alt > 0 && r1Alt < mpr) { mpr = r1Alt; }
          }
        }
        else
        {
          G4double rho, th, rho13, costh3, sinth3, spt, smti32, r1;

          rho = std::sqrt(-q3);
          th = std::acos(r / rho);
          rho13 = std::cbrt(rho);
          costh3 = std::cos(th / 3);
          sinth3 = std::sin(th / 3);
          spt = rho13 * 2 * costh3;
          smti32 = -rho13 * sinth3 * sqrt3;
          r1 = spt - val;
          if (r1 > 0) { mpr = r1; }
          r1 = -spt / 2 - val + smti32;
          if ((r1 > 0) && (r1 < mpr)) { mpr = r1; }
          r1 = r1 - 2 * smti32;
          if ((r1 > 0) && (r1 < mpr)) { mpr = r1; }

          if (discAlt >= 0)
          {
            G4double sdAlt, sAlt, tAlt, r1Alt, rsdAlt;
            sdAlt = std::sqrt(discAlt);
            rsdAlt = rAlt + sdAlt;
            if (rsdAlt > 0) {
              sAlt = std::cbrt(rsdAlt);
            } else {
              sAlt = -std::cbrt(std::fabs(rsdAlt));
            }

            rsdAlt = rAlt - sdAlt;
            if (rsdAlt > 0) {
              tAlt = std::cbrt(rsdAlt);
            } else {
              tAlt = -std::cbrt(std::fabs(rsdAlt));
            }

            r1Alt = sAlt + tAlt - val;

            if (r1Alt > 0 && r1Alt < mpr) { mpr = r1Alt; }
          }
          else
          {
            G4double thAlt, costh3Alt, sinth3Alt, sptAlt, smti32Alt, r1Alt;
            thAlt = std::acos(rAlt / rho);
            costh3Alt = std::cos(thAlt / 3);
            sinth3Alt = std::sin(thAlt / 3);
            sptAlt = rho13 * 2 * costh3Alt;
            smti32Alt = -rho13 * sinth3Alt * sqrt3;
            r1Alt = sptAlt - val;
            if (r1Alt > 0 && r1Alt < mpr) { mpr = r1Alt; }
            r1Alt = -sptAlt / 2 - val + smti32Alt;
            if (r1Alt > 0 && r1Alt < mpr) { mpr = r1Alt; }
            r1Alt = r1Alt - 2 * smti32Alt;
            if (r1Alt > 0 && r1Alt < mpr) { mpr = r1Alt; }
          }
        }
      }

      return mpr;
    }

  private:

    QSS_simulator simulator;
};

#endif
