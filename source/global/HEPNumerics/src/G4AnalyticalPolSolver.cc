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
// G4AnalyticalPolSolver class implementation
//
// Author: V.Grichine, 24.04.97
// --------------------------------------------------------------------

#include "globals.hh"
#include <complex>

#include "G4AnalyticalPolSolver.hh"

//////////////////////////////////////////////////////////////////////////////

G4AnalyticalPolSolver::G4AnalyticalPolSolver() { ; }

//////////////////////////////////////////////////////////////////////////////

G4AnalyticalPolSolver::~G4AnalyticalPolSolver() { ; }

//////////////////////////////////////////////////////////////////////////////
//
// Array r[3][5]  p[5]
// Roots of poly p[0] x^2 + p[1] x+p[2]=0
//
// x = r[1][k] + i r[2][k];  k = 1, 2

G4int G4AnalyticalPolSolver::QuadRoots(G4double p[5], G4double r[3][5])
{
  G4double b, c, d2, d;

  b  = -p[1] / p[0] / 2.;
  c  = p[2] / p[0];
  d2 = b * b - c;

  if(d2 >= 0.)
  {
    d       = std::sqrt(d2);
    r[1][1] = b - d;
    r[1][2] = b + d;
    r[2][1] = 0.;
    r[2][2] = 0.;
  }
  else
  {
    d       = std::sqrt(-d2);
    r[2][1] = d;
    r[2][2] = -d;
    r[1][1] = b;
    r[1][2] = b;
  }

  return 2;
}

//////////////////////////////////////////////////////////////////////////////
//
// Array r[3][5]  p[5]
// Roots of poly p[0] x^3 + p[1] x^2...+p[3]=0
// x=r[1][k] + i r[2][k]  k=1,...,3
// Assumes 0<arctan(x)<pi/2 for x>0

G4int G4AnalyticalPolSolver::CubicRoots(G4double p[5], G4double r[3][5])
{
  G4double x, t, b, c, d;
  G4int k;

  if(p[0] != 1.)
  {
    for(k = 1; k < 4; k++)
    {
      p[k] = p[k] / p[0];
    }
    p[0] = 1.;
  }
  x = p[1] / 3.0;
  t = x * p[1];
  b = 0.5 * (x * (t / 1.5 - p[2]) + p[3]);
  t = (t - p[2]) / 3.0;
  c = t * t * t;
  d = b * b - c;

  if(d >= 0.)
  {
    d = std::pow((std::sqrt(d) + std::fabs(b)), 1.0 / 3.0);

    if(d != 0.)
    {
      if(b > 0.)
      {
        b = -d;
      }
      else
      {
        b = d;
      }
      c = t / b;
    }
    d       = std::sqrt(0.75) * (b - c);
    r[2][2] = d;
    b       = b + c;
    c       = -0.5 * b - x;
    r[1][2] = c;

    if((b > 0. && x <= 0.) || (b < 0. && x > 0.))
    {
      r[1][1] = c;
      r[2][1] = -d;
      r[1][3] = b - x;
      r[2][3] = 0;
    }
    else
    {
      r[1][1] = b - x;
      r[2][1] = 0.;
      r[1][3] = c;
      r[2][3] = -d;
    }
  }  // end of 2 equal or complex roots
  else
  {
    if(b == 0.)
    {
      d = std::atan(1.0) / 1.5;
    }
    else
    {
      d = std::atan(std::sqrt(-d) / std::fabs(b)) / 3.0;
    }

    if(b < 0.)
    {
      b = std::sqrt(t) * 2.0;
    }
    else
    {
      b = -2.0 * std::sqrt(t);
    }

    c = std::cos(d) * b;
    t = -std::sqrt(0.75) * std::sin(d) * b - 0.5 * c;
    d = -t - c - x;
    c = c - x;
    t = t - x;

    if(std::fabs(c) > std::fabs(t))
    {
      r[1][3] = c;
    }
    else
    {
      r[1][3] = t;
      t       = c;
    }
    if(std::fabs(d) > std::fabs(t))
    {
      r[1][2] = d;
    }
    else
    {
      r[1][2] = t;
      t       = d;
    }
    r[1][1] = t;

    for(k = 1; k < 4; k++)
    {
      r[2][k] = 0.;
    }
  }
  return 0;
}

//////////////////////////////////////////////////////////////////////////////
//
// Array r[3][5]  p[5]
// Roots of poly p[0] x^4 + p[1] x^3...+p[4]=0
// x=r[1][k] + i r[2][k]  k=1,...,4

G4int G4AnalyticalPolSolver::BiquadRoots(G4double p[5], G4double r[3][5])
{
  G4double a, b, c, d, e;
  G4int i, k, j;

  if(p[0] != 1.0)
  {
    for(k = 1; k < 5; k++)
    {
      p[k] = p[k] / p[0];
    }
    p[0] = 1.;
  }
  e = 0.25 * p[1];
  b = 2 * e;
  c = b * b;
  d = 0.75 * c;
  b = p[3] + b * (c - p[2]);
  a = p[2] - d;
  c = p[4] + e * (e * a - p[3]);
  a = a - d;

  p[1] = 0.5 * a;
  p[2] = (p[1] * p[1] - c) * 0.25;
  p[3] = b * b / (-64.0);

  if(p[3] < 0.)
  {
    CubicRoots(p, r);

    for(k = 1; k < 4; k++)
    {
      if(r[2][k] == 0. && r[1][k] > 0)
      {
        d = r[1][k] * 4;
        a = a + d;

        if(a >= 0. && b >= 0.)
        {
          p[1] = std::sqrt(d);
        }
        else if(a <= 0. && b <= 0.)
        {
          p[1] = std::sqrt(d);
        }
        else
        {
          p[1] = -std::sqrt(d);
        }

        b = 0.5 * (a + b / p[1]);

        p[2] = c / b;
        QuadRoots(p, r);

        for(i = 1; i < 3; i++)
        {
          for(j = 1; j < 3; j++)
          {
            r[j][i + 2] = r[j][i];
          }
        }
        p[1] = -p[1];
        p[2] = b;
        QuadRoots(p, r);

        for(i = 1; i < 5; i++)
        {
          r[1][i] = r[1][i] - e;
        }

        return 4;
      }
    }
  }
  if(p[2] < 0.)
  {
    b    = std::sqrt(c);
    d    = b + b - a;
    p[1] = 0.;

    if(d > 0.)
    {
      p[1] = std::sqrt(d);
    }
  }
  else
  {
    if(p[1] > 0.)
    {
      b = std::sqrt(p[2]) * 2.0 + p[1];
    }
    else
    {
      b = -std::sqrt(p[2]) * 2.0 + p[1];
    }

    if(b != 0.)
    {
      p[1] = 0;
    }
    else
    {
      for(k = 1; k < 5; k++)
      {
        r[1][k] = -e;
        r[2][k] = 0;
      }
      return 0;
    }
  }

  p[2] = c / b;
  QuadRoots(p, r);

  for(k = 1; k < 3; k++)
  {
    for(j = 1; j < 3; j++)
    {
      r[j][k + 2] = r[j][k];
    }
  }
  p[1] = -p[1];
  p[2] = b;
  QuadRoots(p, r);

  for(k = 1; k < 5; k++)
  {
    r[1][k] = r[1][k] - e;
  }

  return 4;
}

//////////////////////////////////////////////////////////////////////////////

G4int G4AnalyticalPolSolver::QuarticRoots(G4double p[5], G4double r[3][5])
{
  G4double a0, a1, a2, a3, y1;
  G4double R2, D2, E2, D, E, R = 0.;
  G4double a, b, c, d, ds;

  G4double reRoot[4];
  G4int k;

  for(k = 0; k < 4; k++)
  {
    reRoot[k] = DBL_MAX;
  }

  if(p[0] != 1.0)
  {
    for(k = 1; k < 5; k++)
    {
      p[k] = p[k] / p[0];
    }
    p[0] = 1.;
  }
  a3 = p[1];
  a2 = p[2];
  a1 = p[3];
  a0 = p[4];

  // resolvent cubic equation cofs:

  p[1] = -a2;
  p[2] = a1 * a3 - 4 * a0;
  p[3] = 4 * a2 * a0 - a1 * a1 - a3 * a3 * a0;

  CubicRoots(p, r);

  for(k = 1; k < 4; k++)
  {
    if(r[2][k] == 0.)  // find a real root
    {
      reRoot[k] = r[1][k];
    }
    else
    {
      reRoot[k] = DBL_MAX;  // kInfinity;
    }
  }
  y1 = DBL_MAX;  // kInfinity;
  for(k = 1; k < 4; k++)
  {
    if(reRoot[k] < y1)
    {
      y1 = reRoot[k];
    }
  }

  R2 = 0.25 * a3 * a3 - a2 + y1;
  b  = 0.25 * (4 * a3 * a2 - 8 * a1 - a3 * a3 * a3);
  c  = 0.75 * a3 * a3 - 2 * a2;
  a  = c - R2;
  d  = 4 * y1 * y1 - 16 * a0;

  if(R2 > 0.)
  {
    R  = std::sqrt(R2);
    D2 = a + b / R;
    E2 = a - b / R;

    if(D2 >= 0.)
    {
      D       = std::sqrt(D2);
      r[1][1] = -0.25 * a3 + 0.5 * R + 0.5 * D;
      r[1][2] = -0.25 * a3 + 0.5 * R - 0.5 * D;
      r[2][1] = 0.;
      r[2][2] = 0.;
    }
    else
    {
      D       = std::sqrt(-D2);
      r[1][1] = -0.25 * a3 + 0.5 * R;
      r[1][2] = -0.25 * a3 + 0.5 * R;
      r[2][1] = 0.5 * D;
      r[2][2] = -0.5 * D;
    }
    if(E2 >= 0.)
    {
      E       = std::sqrt(E2);
      r[1][3] = -0.25 * a3 - 0.5 * R + 0.5 * E;
      r[1][4] = -0.25 * a3 - 0.5 * R - 0.5 * E;
      r[2][3] = 0.;
      r[2][4] = 0.;
    }
    else
    {
      E       = std::sqrt(-E2);
      r[1][3] = -0.25 * a3 - 0.5 * R;
      r[1][4] = -0.25 * a3 - 0.5 * R;
      r[2][3] = 0.5 * E;
      r[2][4] = -0.5 * E;
    }
  }
  else if(R2 < 0.)
  {
    R = std::sqrt(-R2);
    G4complex CD2(a, -b / R);
    G4complex CD = std::sqrt(CD2);

    r[1][1] = -0.25 * a3 + 0.5 * real(CD);
    r[1][2] = -0.25 * a3 - 0.5 * real(CD);
    r[2][1] = 0.5 * R + 0.5 * imag(CD);
    r[2][2] = 0.5 * R - 0.5 * imag(CD);
    G4complex CE2(a, b / R);
    G4complex CE = std::sqrt(CE2);

    r[1][3] = -0.25 * a3 + 0.5 * real(CE);
    r[1][4] = -0.25 * a3 - 0.5 * real(CE);
    r[2][3] = -0.5 * R + 0.5 * imag(CE);
    r[2][4] = -0.5 * R - 0.5 * imag(CE);
  }
  else  // R2=0 case
  {
    if(d >= 0.)
    {
      D2 = c + std::sqrt(d);
      E2 = c - std::sqrt(d);

      if(D2 >= 0.)
      {
        D       = std::sqrt(D2);
        r[1][1] = -0.25 * a3 + 0.5 * R + 0.5 * D;
        r[1][2] = -0.25 * a3 + 0.5 * R - 0.5 * D;
        r[2][1] = 0.;
        r[2][2] = 0.;
      }
      else
      {
        D       = std::sqrt(-D2);
        r[1][1] = -0.25 * a3 + 0.5 * R;
        r[1][2] = -0.25 * a3 + 0.5 * R;
        r[2][1] = 0.5 * D;
        r[2][2] = -0.5 * D;
      }
      if(E2 >= 0.)
      {
        E       = std::sqrt(E2);
        r[1][3] = -0.25 * a3 - 0.5 * R + 0.5 * E;
        r[1][4] = -0.25 * a3 - 0.5 * R - 0.5 * E;
        r[2][3] = 0.;
        r[2][4] = 0.;
      }
      else
      {
        E       = std::sqrt(-E2);
        r[1][3] = -0.25 * a3 - 0.5 * R;
        r[1][4] = -0.25 * a3 - 0.5 * R;
        r[2][3] = 0.5 * E;
        r[2][4] = -0.5 * E;
      }
    }
    else
    {
      ds = std::sqrt(-d);
      G4complex CD2(c, ds);
      G4complex CD = std::sqrt(CD2);

      r[1][1] = -0.25 * a3 + 0.5 * real(CD);
      r[1][2] = -0.25 * a3 - 0.5 * real(CD);
      r[2][1] = 0.5 * R + 0.5 * imag(CD);
      r[2][2] = 0.5 * R - 0.5 * imag(CD);

      G4complex CE2(c, -ds);
      G4complex CE = std::sqrt(CE2);

      r[1][3] = -0.25 * a3 + 0.5 * real(CE);
      r[1][4] = -0.25 * a3 - 0.5 * real(CE);
      r[2][3] = -0.5 * R + 0.5 * imag(CE);
      r[2][4] = -0.5 * R - 0.5 * imag(CE);
    }
  }
  return 4;
}

//
//
//////////////////////////////////////////////////////////////////////////////
