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
// G4Exp
//
// Class description:
//
// The basic idea is to exploit Pade polynomials.
// A lot of ideas were inspired by the cephes math library
// (by Stephen L. Moshier moshier@na-net.ornl.gov) as well as actual code.
// The Cephes library can be found here:  http://www.netlib.org/cephes/
// Code and algorithms for G4Exp have been extracted and adapted for Geant4
// from the original implementation in the VDT mathematical library
// (https://svnweb.cern.ch/trac/vdt), version 0.3.7.

// Original implementation created on: Jun 23, 2012
// Authors: Danilo Piparo, Thomas Hauth, Vincenzo Innocente
//
// --------------------------------------------------------------------
/*
 * VDT is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser Public License for more details.
 *
 * You should have received a copy of the GNU Lesser Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */
// --------------------------------------------------------------------
#ifndef G4Exp_hh
#define G4Exp_hh 1

#ifdef WIN32

#  define G4Exp std::exp

#else

#  include "G4Types.hh"

#  include <cstdint>
#  include <limits>

namespace G4ExpConsts
{
  const G4double EXP_LIMIT = 708;

  const G4double PX1exp = 1.26177193074810590878E-4;
  const G4double PX2exp = 3.02994407707441961300E-2;
  const G4double PX3exp = 9.99999999999999999910E-1;
  const G4double QX1exp = 3.00198505138664455042E-6;
  const G4double QX2exp = 2.52448340349684104192E-3;
  const G4double QX3exp = 2.27265548208155028766E-1;
  const G4double QX4exp = 2.00000000000000000009E0;

  const G4double LOG2E = 1.4426950408889634073599;  // 1/log(2)

  const G4float MAXLOGF = 88.72283905206835f;
  const G4float MINLOGF = -88.f;

  const G4float C1F = 0.693359375f;
  const G4float C2F = -2.12194440e-4f;

  const G4float PX1expf = 1.9875691500E-4f;
  const G4float PX2expf = 1.3981999507E-3f;
  const G4float PX3expf = 8.3334519073E-3f;
  const G4float PX4expf = 4.1665795894E-2f;
  const G4float PX5expf = 1.6666665459E-1f;
  const G4float PX6expf = 5.0000001201E-1f;

  const G4float LOG2EF = 1.44269504088896341f;

  //----------------------------------------------------------------------------
  // Used to switch between different type of interpretations of the data
  // (64 bits)
  //
  union ieee754
  {
    ieee754()= default;
    ieee754(G4double thed) { d = thed; };
    ieee754(uint64_t thell) { ll = thell; };
    ieee754(G4float thef) { f[0] = thef; };
    ieee754(uint32_t thei) { i[0] = thei; };
    G4double d;
    G4float f[2];
    uint32_t i[2];
    uint64_t ll;
    uint16_t s[4];
  };

  //----------------------------------------------------------------------------
  // Converts an unsigned long long to a double
  //
  inline G4double uint642dp(uint64_t ll)
  {
    ieee754 tmp;
    tmp.ll = ll;
    return tmp.d;
  }

  //----------------------------------------------------------------------------
  // Converts an int to a float
  //
  inline G4float uint322sp(G4int x)
  {
    ieee754 tmp;
    tmp.i[0] = x;
    return tmp.f[0];
  }

  //----------------------------------------------------------------------------
  // Converts a float to an int
  //
  inline uint32_t sp2uint32(G4float x)
  {
    ieee754 tmp;
    tmp.f[0] = x;
    return tmp.i[0];
  }

  //----------------------------------------------------------------------------
  /**
   * A vectorisable floor implementation, not only triggered by fast-math.
   * These functions do not distinguish between -0.0 and 0.0, so are not IEC6509
   * compliant for argument -0.0
   **/
  inline G4double fpfloor(const G4double x)
  {
    // no problem since exp is defined between -708 and 708. Int is enough for
    // it!
    int32_t ret = int32_t(x);
    ret -= (sp2uint32(x) >> 31);
    return ret;
  }

  //----------------------------------------------------------------------------
  /**
   * A vectorisable floor implementation, not only triggered by fast-math.
   * These functions do not distinguish between -0.0 and 0.0, so are not IEC6509
   * compliant for argument -0.0
   **/
  inline G4float fpfloor(const G4float x)
  {
    int32_t ret = int32_t(x);
    ret -= (sp2uint32(x) >> 31);
    return ret;
  }
}  // namespace G4ExpConsts

// Exp double precision --------------------------------------------------------

/// Exponential Function double precision
inline G4double G4Exp(G4double initial_x)
{
  G4double x  = initial_x;
  G4double px = G4ExpConsts::fpfloor(G4ExpConsts::LOG2E * x + 0.5);

  const int32_t n = int32_t(px);

  x -= px * 6.93145751953125E-1;
  x -= px * 1.42860682030941723212E-6;

  const G4double xx = x * x;

  // px = x * P(x**2).
  px = G4ExpConsts::PX1exp;
  px *= xx;
  px += G4ExpConsts::PX2exp;
  px *= xx;
  px += G4ExpConsts::PX3exp;
  px *= x;

  // Evaluate Q(x**2).
  G4double qx = G4ExpConsts::QX1exp;
  qx *= xx;
  qx += G4ExpConsts::QX2exp;
  qx *= xx;
  qx += G4ExpConsts::QX3exp;
  qx *= xx;
  qx += G4ExpConsts::QX4exp;

  // e**x = 1 + 2x P(x**2)/( Q(x**2) - P(x**2) )
  x = px / (qx - px);
  x = 1.0 + 2.0 * x;

  // Build 2^n in double.
  x *= G4ExpConsts::uint642dp((((uint64_t) n) + 1023) << 52);

  if(initial_x > G4ExpConsts::EXP_LIMIT)
    x = std::numeric_limits<G4double>::infinity();
  if(initial_x < -G4ExpConsts::EXP_LIMIT)
    x = 0.;

  return x;
}

// Exp single precision --------------------------------------------------------

/// Exponential Function single precision
inline G4float G4Expf(G4float initial_x)
{
  G4float x = initial_x;

  G4float z =
    G4ExpConsts::fpfloor(G4ExpConsts::LOG2EF * x +
                         0.5f); /* std::floor() truncates toward -infinity. */

  x -= z * G4ExpConsts::C1F;
  x -= z * G4ExpConsts::C2F;
  const int32_t n = int32_t(z);

  const G4float x2 = x * x;

  z = x * G4ExpConsts::PX1expf;
  z += G4ExpConsts::PX2expf;
  z *= x;
  z += G4ExpConsts::PX3expf;
  z *= x;
  z += G4ExpConsts::PX4expf;
  z *= x;
  z += G4ExpConsts::PX5expf;
  z *= x;
  z += G4ExpConsts::PX6expf;
  z *= x2;
  z += x + 1.0f;

  /* multiply by power of 2 */
  z *= G4ExpConsts::uint322sp((n + 0x7f) << 23);

  if(initial_x > G4ExpConsts::MAXLOGF)
    z = std::numeric_limits<G4float>::infinity();
  if(initial_x < G4ExpConsts::MINLOGF)
    z = 0.f;

  return z;
}

//------------------------------------------------------------------------------

void expv(const uint32_t size, G4double const* __restrict__ iarray,
          G4double* __restrict__ oarray);
void G4Expv(const uint32_t size, G4double const* __restrict__ iarray,
            G4double* __restrict__ oarray);
void expfv(const uint32_t size, G4float const* __restrict__ iarray,
           G4float* __restrict__ oarray);
void G4Expfv(const uint32_t size, G4float const* __restrict__ iarray,
             G4float* __restrict__ oarray);

#endif /* WIN32 */

#endif
