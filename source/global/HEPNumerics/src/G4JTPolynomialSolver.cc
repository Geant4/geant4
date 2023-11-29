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
// G4JTPolynomialSolver class implementation.
//
// Author: Oliver Link, 15.02.2005
//         Translated to C++ and adapted to use STL vectors.
// --------------------------------------------------------------------

#include "G4JTPolynomialSolver.hh"
#include "G4Pow.hh"
#include "G4SystemOfUnits.hh"

const G4double G4JTPolynomialSolver::base   = 2;
const G4double G4JTPolynomialSolver::eta    = DBL_EPSILON;
const G4double G4JTPolynomialSolver::infin  = DBL_MAX;
const G4double G4JTPolynomialSolver::smalno = DBL_MIN;
const G4double G4JTPolynomialSolver::are    = DBL_EPSILON;
const G4double G4JTPolynomialSolver::mre    = DBL_EPSILON;
const G4double G4JTPolynomialSolver::lo     = DBL_MIN / DBL_EPSILON;

G4int G4JTPolynomialSolver::FindRoots(G4double* op, G4int degr, G4double* zeror,
                                      G4double* zeroi)
{
  G4double t = 0.0, aa = 0.0, bb = 0.0, cc = 0.0, factor = 1.0;
  G4double max = 0.0, min = infin, xxx = 0.0, x = 0.0, sc = 0.0, bnd = 0.0;
  G4double xm = 0.0, ff = 0.0, df = 0.0, dx = 0.0;
  G4int cnt = 0, nz = 0, i = 0, j = 0, jj = 0, l = 0, nm1 = 0, zerok = 0;
  G4Pow* power = G4Pow::GetInstance();

  // Initialization of constants for shift rotation.
  //
  static const G4double xx   = std::sqrt(0.5);
  static const G4double rot  = 94.0 * deg;
  static const G4double cosr = std::cos(rot), sinr = std::sin(rot);
  G4double xo = xx, yo = -xx;
  n = degr;

  //  Algorithm fails if the leading coefficient is zero.
  //
  if(!(op[0] != 0.0))
  {
    return -1;
  }

  //  Remove the zeros at the origin, if any.
  //
  while(!(op[n] != 0.0))
  {
    j        = degr - n;
    zeror[j] = 0.0;
    zeroi[j] = 0.0;
    n--;
  }
  if(n < 1)
  {
    return -1;
  }

  // Allocate buffers here
  //
  std::vector<G4double> temp(degr + 1);
  std::vector<G4double> pt(degr + 1);

  p.assign(degr + 1, 0);
  qp.assign(degr + 1, 0);
  k.assign(degr + 1, 0);
  qk.assign(degr + 1, 0);
  svk.assign(degr + 1, 0);

  // Make a copy of the coefficients.
  //
  for(i = 0; i <= n; ++i)
  {
    p[i] = op[i];
  }

  do
  {
    if(n == 1)  // Start the algorithm for one zero.
    {
      zeror[degr - 1] = -p[1] / p[0];
      zeroi[degr - 1] = 0.0;
      n -= 1;
      return degr - n;
    }
    if(n == 2)  // Calculate the final zero or pair of zeros.
    {
      Quadratic(p[0], p[1], p[2], &zeror[degr - 2], &zeroi[degr - 2],
                &zeror[degr - 1], &zeroi[degr - 1]);
      n -= 2;
      return degr - n;
    }

    // Find largest and smallest moduli of coefficients.
    //
    max = 0.0;
    min = infin;
    for(i = 0; i <= n; ++i)
    {
      x = std::fabs(p[i]);
      if(x > max)
      {
        max = x;
      }
      if(x != 0.0 && x < min)
      {
        min = x;
      }
    }

    // Scale if there are large or very small coefficients.
    // Computes a scale factor to multiply the coefficients of the
    // polynomial. The scaling is done to avoid overflow and to
    // avoid undetected underflow interfering with the convergence
    // criterion. The factor is a power of the base.
    //
    sc = lo / min;

    if(((sc <= 1.0) && (max >= 10.0)) || ((sc > 1.0) && (infin / sc >= max)) ||
       ((infin / sc >= max) && (max >= 10)))
    {
      if(!(sc != 0.0))
      {
        sc = smalno;
      }
      l      = (G4int)(G4Log(sc) / G4Log(base) + 0.5);
      factor = power->powN(base, l);
      if(factor != 1.0)
      {
        for(i = 0; i <= n; ++i)
        {
          p[i] = factor * p[i];
        }  // Scale polynomial.
      }
    }

    // Compute lower bound on moduli of roots.
    //
    for(i = 0; i <= n; ++i)
    {
      pt[i] = (std::fabs(p[i]));
    }
    pt[n] = -pt[n];

    // Compute upper estimate of bound.
    //
    x = G4Exp((G4Log(-pt[n]) - G4Log(pt[0])) / (G4double) n);

    // If Newton step at the origin is better, use it.
    //
    if(pt[n - 1] != 0.0)
    {
      xm = -pt[n] / pt[n - 1];
      if(xm < x)
      {
        x = xm;
      }
    }

    // Chop the interval (0,x) until ff <= 0
    //
    while(true)
    {
      xm = x * 0.1;
      ff = pt[0];
      for(i = 1; i <= n; ++i)
      {
        ff = ff * xm + pt[i];
      }
      if(ff <= 0.0)
      {
        break;
      }
      x = xm;
    }
    dx = x;

    // Do Newton interation until x converges to two decimal places.
    //
    while(std::fabs(dx / x) > 0.005)
    {
      ff = pt[0];
      df = ff;
      for(i = 1; i < n; ++i)
      {
        ff = ff * x + pt[i];
        df = df * x + ff;
      }
      ff = ff * x + pt[n];
      dx = ff / df;
      x -= dx;
    }
    bnd = x;

    // Compute the derivative as the initial k polynomial
    // and do 5 steps with no shift.
    //
    nm1 = n - 1;
    for(i = 1; i < n; ++i)
    {
      k[i] = (G4double)(n - i) * p[i] / (G4double) n;
    }
    k[0]  = p[0];
    aa    = p[n];
    bb    = p[n - 1];
    zerok = static_cast<G4int>(k[n - 1] == 0);
    for(jj = 0; jj < 5; ++jj)
    {
      cc = k[n - 1];
      if(zerok == 0)  // Use a scaled form of recurrence if k at 0 is nonzero.
      {
        // Use a scaled form of recurrence if value of k at 0 is nonzero.
        //
        t = -aa / cc;
        for(i = 0; i < nm1; ++i)
        {
          j    = n - i - 1;
          k[j] = t * k[j - 1] + p[j];
        }
        k[0]  = p[0];
        zerok =
          static_cast<G4int>(std::fabs(k[n - 1]) <= std::fabs(bb) * eta * 10.0);
      }
      else  // Use unscaled form of recurrence.
      {
        for(i = 0; i < nm1; ++i)
        {
          j    = n - i - 1;
          k[j] = k[j - 1];
        }
        k[0]  = 0.0;
        zerok = static_cast<G4int>(!(k[n - 1] != 0.0));
      }
    }

    // Save k for restarts with new shifts.
    //
    for(i = 0; i < n; ++i)
    {
      temp[i] = k[i];
    }

    // Loop to select the quadratic corresponding to each new shift.
    //
    for(cnt = 0; cnt < 20; ++cnt)
    {
      // Quadratic corresponds to a double shift to a
      // non-real point and its complex conjugate. The point
      // has modulus bnd and amplitude rotated by 94 degrees
      // from the previous shift.
      //
      xxx = cosr * xo - sinr * yo;
      yo  = sinr * xo + cosr * yo;
      xo  = xxx;
      sr  = bnd * xo;
      si  = bnd * yo;
      u   = -2.0 * sr;
      v   = bnd;
      ComputeFixedShiftPolynomial(20 * (cnt + 1), &nz);
      if(nz != 0)
      {
        // The second stage jumps directly to one of the third
        // stage iterations and returns here if successful.
        // Deflate the polynomial, store the zero or zeros and
        // return to the main algorithm.
        //
        j        = degr - n;
        zeror[j] = szr;
        zeroi[j] = szi;
        n -= nz;
        for(i = 0; i <= n; ++i)
        {
          p[i] = qp[i];
        }
        if(nz != 1)
        {
          zeror[j + 1] = lzr;
          zeroi[j + 1] = lzi;
        }
        break;
      }
      
      // If the iteration is unsuccessful another quadratic
      // is chosen after restoring k.
      //
      for(i = 0; i < n; ++i)
      {
        k[i] = temp[i];
      }
     
    }
  } while(nz != 0);  // End of initial DO loop

  // Return with failure if no convergence with 20 shifts.
  //
  return degr - n;
}

void G4JTPolynomialSolver::ComputeFixedShiftPolynomial(G4int l2, G4int* nz)
{
  // Computes up to L2 fixed shift k-polynomials, testing for convergence
  // in the linear or quadratic case. Initiates one of the variable shift
  // iterations and returns with the number of zeros found.

  G4double svu = 0.0, svv = 0.0, ui = 0.0, vi = 0.0, xs = 0.0;
  G4double betas = 0.25, betav = 0.25, oss = sr, ovv = v, ss = 0.0, vv = 0.0,
           ts = 1.0, tv = 1.0;
  G4double ots = 0.0, otv = 0.0;
  G4double tvv = 1.0, tss = 1.0;
  G4int type = 0, i = 0, j = 0, iflag = 0, vpass = 0, spass = 0, vtry = 0,
        stry = 0;

  *nz = 0;

  // Evaluate polynomial by synthetic division.
  //
  QuadraticSyntheticDivision(n, &u, &v, p, qp, &a, &b);
  ComputeScalarFactors(&type);
  for(j = 0; j < l2; ++j)
  {
    // Calculate next k polynomial and estimate v.
    //
    ComputeNextPolynomial(&type);
    ComputeScalarFactors(&type);
    ComputeNewEstimate(type, &ui, &vi);
    vv = vi;

    // Estimate xs.
    //
    ss = 0.0;
    if(k[n - 1] != 0.0)
    {
      ss = -p[n] / k[n - 1];
    }
    tv = 1.0;
    ts = 1.0;
    if(j == 0 || type == 3)
    {
      ovv = vv;
      oss = ss;
      otv = tv;
      ots = ts;
      continue;
    }

    // Compute relative measures of convergence of xs and v sequences.
    //
    if(vv != 0.0)
    {
      tv = std::fabs((vv - ovv) / vv);
    }
    if(ss != 0.0)
    {
      ts = std::fabs((ss - oss) / ss);
    }

    // If decreasing, multiply two most recent convergence measures.
    tvv = 1.0;
    if(tv < otv)
    {
      tvv = tv * otv;
    }
    tss = 1.0;
    if(ts < ots)
    {
      tss = ts * ots;
    }

    // Compare with convergence criteria.
    vpass = static_cast<G4int>(tvv < betav);
    spass = static_cast<G4int>(tss < betas);
    if(!((spass != 0) || (vpass != 0)))
    {
      ovv = vv;
      oss = ss;
      otv = tv;
      ots = ts;
      continue;
    }

    // At least one sequence has passed the convergence test.
    // Store variables before iterating.
    //
    svu = u;
    svv = v;
    for(i = 0; i < n; ++i)
    {
      svk[i] = k[i];
    }
    xs = ss;

    // Choose iteration according to the fastest converging sequence.
    //
    vtry = 0;
    stry = 0;
    if(((spass != 0) && (vpass == 0)) || (tss < tvv))
    {
      RealPolynomialIteration(&xs, nz, &iflag);
      if(*nz > 0)
      {
        return;
      }

      // Linear iteration has failed. Flag that it has been
      // tried and decrease the convergence criterion.
      //
      stry = 1;
      betas *= 0.25;
      if(iflag == 0)
      {
        goto _restore_variables;
      }

      // If linear iteration signals an almost double real
      // zero attempt quadratic iteration.
      //
      ui = -(xs + xs);
      vi = xs * xs;
    }

  _quadratic_iteration:

    do
    {
      QuadraticPolynomialIteration(&ui, &vi, nz);
      if(*nz > 0)
      {
        return;
      }

      // Quadratic iteration has failed. Flag that it has
      // been tried and decrease the convergence criterion.
      //
      vtry = 1;
      betav *= 0.25;

      // Try linear iteration if it has not been tried and
      // the S sequence is converging.
      //
      if((stry != 0) || (spass == 0))
      {
        break;
      }
      for(i = 0; i < n; ++i)
      {
        k[i] = svk[i];
      }
      RealPolynomialIteration(&xs, nz, &iflag);
      if(*nz > 0)
      {
        return;
      }

      // Linear iteration has failed. Flag that it has been
      // tried and decrease the convergence criterion.
      //
      stry = 1;
      betas *= 0.25;
      if(iflag == 0)
      {
        break;
      }

      // If linear iteration signals an almost double real
      // zero attempt quadratic iteration.
      //
      ui = -(xs + xs);
      vi = xs * xs;
    } while(iflag != 0);

    // Restore variables.

  _restore_variables:

    u = svu;
    v = svv;
    for(i = 0; i < n; ++i)
    {
      k[i] = svk[i];
    }

    // Try quadratic iteration if it has not been tried
    // and the V sequence is converging.
    //
    if((vpass != 0) && (vtry == 0))
    {
      goto _quadratic_iteration;
    }

    // Recompute QP and scalar values to continue the
    // second stage.
    //
    QuadraticSyntheticDivision(n, &u, &v, p, qp, &a, &b);
    ComputeScalarFactors(&type);

    ovv = vv;
    oss = ss;
    otv = tv;
    ots = ts;
  }
}

void G4JTPolynomialSolver::QuadraticPolynomialIteration(G4double* uu,
                                                        G4double* vv, G4int* nz)
{
  // Variable-shift k-polynomial iteration for a
  // quadratic factor converges only if the zeros are
  // equimodular or nearly so.
  // uu, vv - coefficients of starting quadratic.
  // nz - number of zeros found.
  //
  G4double ui = 0.0, vi = 0.0;
  G4double omp    = 0.0;
  G4double relstp = 0.0;
  G4double mp = 0.0, ee = 0.0, t = 0.0, zm = 0.0;
  G4int type = 0, i = 1, j = 0, tried = 0;

  *nz   = 0;
  tried = 0;
  u     = *uu;
  v     = *vv;

  // Main loop.

  while(true)
  {
    Quadratic(1.0, u, v, &szr, &szi, &lzr, &lzi);

    // Return if roots of the quadratic are real and not
    // close to multiple or nearly equal and of opposite
    // sign.
    //
    if(std::fabs(std::fabs(szr) - std::fabs(lzr)) > 0.01 * std::fabs(lzr))
    {
      return;
    }

    // Evaluate polynomial by quadratic synthetic division.
    //
    QuadraticSyntheticDivision(n, &u, &v, p, qp, &a, &b);
    mp = std::fabs(a - szr * b) + std::fabs(szi * b);

    // Compute a rigorous bound on the rounding error in evaluating p.
    //
    zm = std::sqrt(std::fabs(v));
    ee = 2.0 * std::fabs(qp[0]);
    t  = -szr * b;
    for(i = 1; i < n; ++i)
    {
      ee = ee * zm + std::fabs(qp[i]);
    }
    ee = ee * zm + std::fabs(a + t);
    ee *= (5.0 * mre + 4.0 * are);
    ee = ee - (5.0 * mre + 2.0 * are) * (std::fabs(a + t) + std::fabs(b) * zm) +
         2.0 * are * std::fabs(t);

    // Iteration has converged sufficiently if the
    // polynomial value is less than 20 times this bound.
    //
    if(mp <= 20.0 * ee)
    {
      *nz = 2;
      return;
    }
    j++;

    // Stop iteration after 20 steps.
    //
    if(j > 20)
    {
      return;
    }
    if(j >= 2)
    {
      if(!(relstp > 0.01 || mp < omp || (tried != 0)))
      {
        // A cluster appears to be stalling the convergence.
        // Five fixed shift steps are taken with a u,v close to the cluster.
        //
        if(relstp < eta)
        {
          relstp = eta;
        }
        relstp = std::sqrt(relstp);
        u      = u - u * relstp;
        v      = v + v * relstp;
        QuadraticSyntheticDivision(n, &u, &v, p, qp, &a, &b);
        for(i = 0; i < 5; ++i)
        {
          ComputeScalarFactors(&type);
          ComputeNextPolynomial(&type);
        }
        tried = 1;
        j     = 0;
      }
    }
    omp = mp;

    // Calculate next k polynomial and new u and v.
    //
    ComputeScalarFactors(&type);
    ComputeNextPolynomial(&type);
    ComputeScalarFactors(&type);
    ComputeNewEstimate(type, &ui, &vi);

    // If vi is zero the iteration is not converging.
    //
    if(!(vi != 0.0))
    {
      return;
    }
    relstp = std::fabs((vi - v) / vi);
    u      = ui;
    v      = vi;
  }
}

void G4JTPolynomialSolver::RealPolynomialIteration(G4double* sss, G4int* nz,
                                                   G4int* iflag)
{
  // Variable-shift H polynomial iteration for a real zero.
  // sss - starting iterate
  // nz  - number of zeros found
  // iflag - flag to indicate a pair of zeros near real axis.

  G4double t   = 0.;
  G4double omp = 0.;
  G4double pv = 0.0, kv = 0.0, xs = *sss;
  G4double mx = 0.0, mp = 0.0, ee = 0.0;
  G4int i = 1, j = 0;

  *nz    = 0;
  *iflag = 0;

  // Main loop
  //
  while(true)
  {
    pv = p[0];

    // Evaluate p at xs.
    //
    qp[0] = pv;
    for(i = 1; i <= n; ++i)
    {
      pv    = pv * xs + p[i];
      qp[i] = pv;
    }
    mp = std::fabs(pv);

    // Compute a rigorous bound on the error in evaluating p.
    //
    mx = std::fabs(xs);
    ee = (mre / (are + mre)) * std::fabs(qp[0]);
    for(i = 1; i <= n; ++i)
    {
      ee = ee * mx + std::fabs(qp[i]);
    }

    // Iteration has converged sufficiently if the polynomial
    // value is less than 20 times this bound.
    //
    if(mp <= 20.0 * ((are + mre) * ee - mre * mp))
    {
      *nz = 1;
      szr = xs;
      szi = 0.0;
      return;
    }
    j++;

    // Stop iteration after 10 steps.
    //
    if(j > 10)
    {
      return;
    }
    if(j >= 2)
    {
      if(!(std::fabs(t) > 0.001 * std::fabs(xs - t) || mp < omp))
      {
        // A cluster of zeros near the real axis has been encountered.
        // Return with iflag set to initiate a quadratic iteration.
        //
        *iflag = 1;
        *sss   = xs;
        return;
      }  // Return if the polynomial value has increased significantly.
    }

    omp = mp;

    //  Compute t, the next polynomial, and the new iterate.
    //
    kv    = k[0];
    qk[0] = kv;
    for(i = 1; i < n; ++i)
    {
      kv    = kv * xs + k[i];
      qk[i] = kv;
    }
    if(std::fabs(kv) <= std::fabs(k[n - 1]) * 10.0 * eta)  // Use unscaled form.
    {
      k[0] = 0.0;
      for(i = 1; i < n; ++i)
      {
        k[i] = qk[i - 1];
      }
    }
    else  // Use the scaled form of the recurrence if k at xs is nonzero.
    {
      t    = -pv / kv;
      k[0] = qp[0];
      for(i = 1; i < n; ++i)
      {
        k[i] = t * qk[i - 1] + qp[i];
      }
    }
    kv = k[0];
    for(i = 1; i < n; ++i)
    {
      kv = kv * xs + k[i];
    }
    t = 0.0;
    if(std::fabs(kv) > std::fabs(k[n - 1] * 10.0 * eta))
    {
      t = -pv / kv;
    }
    xs += t;
  }
}

void G4JTPolynomialSolver::ComputeScalarFactors(G4int* type)
{
  // This function calculates scalar quantities used to
  // compute the next k polynomial and new estimates of
  // the quadratic coefficients.
  // type - integer variable set here indicating how the
  // calculations are normalized to avoid overflow.

  //  Synthetic division of k by the quadratic 1,u,v
  //
  QuadraticSyntheticDivision(n - 1, &u, &v, k, qk, &c, &d);
  if(std::fabs(c) <= std::fabs(k[n - 1] * 100.0 * eta))
  {
    if(std::fabs(d) <= std::fabs(k[n - 2] * 100.0 * eta))
    {
      *type = 3;  // Type=3 indicates the quadratic is almost a factor of k.
      return;
    }
  }

  if(std::fabs(d) < std::fabs(c))
  {
    *type = 1;  // Type=1 indicates that all formulas are divided by c.
    e     = a / c;
    f     = d / c;
    g     = u * e;
    h     = v * b;
    a3    = a * e + (h / c + g) * b;
    a1    = b - a * (d / c);
    a7    = a + g * d + h * f;
    return;
  }
  *type = 2;  // Type=2 indicates that all formulas are divided by d.
  e     = a / d;
  f     = c / d;
  g     = u * b;
  h     = v * b;
  a3    = (a + g) * e + h * (b / d);
  a1    = b * f - a;
  a7    = (f + u) * a + h;
}

void G4JTPolynomialSolver::ComputeNextPolynomial(G4int* type)
{
  // Computes the next k polynomials using scalars
  // computed in ComputeScalarFactors.

  G4int i = 2;

  if(*type == 3)  // Use unscaled form of the recurrence if type is 3.
  {
    k[0] = 0.0;
    k[1] = 0.0;
    for(i = 2; i < n; ++i)
    {
      k[i] = qk[i - 2];
    }
    return;
  }
  G4double temp = a;
  if(*type == 1)
  {
    temp = b;
  }
  if(std::fabs(a1) <= std::fabs(temp) * eta * 10.0)
  {
    // If a1 is nearly zero then use a special form of the recurrence.
    //
    k[0] = 0.0;
    k[1] = -a7 * qp[0];
    for(i = 2; i < n; ++i)
    {
      k[i] = a3 * qk[i - 2] - a7 * qp[i - 1];
    }
    return;
  }

  // Use scaled form of the recurrence.
  //
  a7 /= a1;
  a3 /= a1;
  k[0] = qp[0];
  k[1] = qp[1] - a7 * qp[0];
  for(i = 2; i < n; ++i)
  {
    k[i] = a3 * qk[i - 2] - a7 * qp[i - 1] + qp[i];
  }
}

void G4JTPolynomialSolver::ComputeNewEstimate(G4int type, G4double* uu,
                                              G4double* vv)
{
  // Compute new estimates of the quadratic coefficients
  // using the scalars computed in calcsc.

  G4double a4 = 0.0, a5 = 0.0, b1 = 0.0, b2 = 0.0, c1 = 0.0, c2 = 0.0, c3 = 0.0,
           c4 = 0.0, temp = 0.0;

  // Use formulas appropriate to setting of type.
  //
  if(type == 3)  //  If type=3 the quadratic is zeroed.
  {
    *uu = 0.0;
    *vv = 0.0;
    return;
  }
  if(type == 2)
  {
    a4 = (a + g) * f + h;
    a5 = (f + u) * c + v * d;
  }
  else
  {
    a4 = a + u * b + h * f;
    a5 = c + (u + v * f) * d;
  }

  //  Evaluate new quadratic coefficients.
  //
  b1   = -k[n - 1] / p[n];
  b2   = -(k[n - 2] + b1 * p[n - 1]) / p[n];
  c1   = v * b2 * a1;
  c2   = b1 * a7;
  c3   = b1 * b1 * a3;
  c4   = c1 - c2 - c3;
  temp = a5 + b1 * a4 - c4;
  if(!(temp != 0.0))
  {
    *uu = 0.0;
    *vv = 0.0;
    return;
  }
  *uu = u - (u * (c3 + c2) + v * (b1 * a1 + b2 * a7)) / temp;
  *vv = v * (1.0 + c4 / temp);
  return;
}

void G4JTPolynomialSolver::QuadraticSyntheticDivision(
  G4int nn, G4double* uu, G4double* vv, std::vector<G4double>& pp,
  std::vector<G4double>& qq, G4double* aa, G4double* bb)
{
  // Divides pp by the quadratic 1,uu,vv placing the quotient
  // in qq and the remainder in aa,bb.

  G4double cc = 0.0;
  *bb         = pp[0];
  qq[0]       = *bb;
  *aa         = pp[1] - (*bb) * (*uu);
  qq[1]       = *aa;
  for(G4int i = 2; i <= nn; ++i)
  {
    cc    = pp[i] - (*aa) * (*uu) - (*bb) * (*vv);
    qq[i] = cc;
    *bb   = *aa;
    *aa   = cc;
  }
}

void G4JTPolynomialSolver::Quadratic(G4double aa, G4double b1, G4double cc,
                                     G4double* ssr, G4double* ssi, G4double* lr,
                                     G4double* li)
{
  // Calculate the zeros of the quadratic aa*z^2 + b1*z + cc.
  // The quadratic formula, modified to avoid overflow, is used
  // to find the larger zero if the zeros are real and both
  // are complex. The smaller real zero is found directly from
  // the product of the zeros c/a.

  G4double bb = 0.0, dd = 0.0, ee = 0.0;

  if(!(aa != 0.0))  // less than two roots
  {
    if(b1 != 0.0)
    {
      *ssr = -cc / b1;
    }
    else
    {
      *ssr = 0.0;
    }
    *lr  = 0.0;
    *ssi = 0.0;
    *li  = 0.0;
    return;
  }
  if(!(cc != 0.0))  // one real root, one zero root
  {
    *ssr = 0.0;
    *lr  = -b1 / aa;
    *ssi = 0.0;
    *li  = 0.0;
    return;
  }

  // Compute discriminant avoiding overflow.
  //
  bb = b1 / 2.0;
  if(std::fabs(bb) < std::fabs(cc))
  {
    if(cc < 0.0)
    {
      ee = -aa;
    }
    else
    {
      ee = aa;
    }
    ee = bb * (bb / std::fabs(cc)) - ee;
    dd = std::sqrt(std::fabs(ee)) * std::sqrt(std::fabs(cc));
  }
  else
  {
    ee = 1.0 - (aa / bb) * (cc / bb);
    dd = std::sqrt(std::fabs(ee)) * std::fabs(bb);
  }
  if(ee < 0.0)  // complex conjugate zeros
  {
    *ssr = -bb / aa;
    *lr  = *ssr;
    *ssi = std::fabs(dd / aa);
    *li  = -(*ssi);
  }
  else
  {
    if(bb >= 0.0)  // real zeros.
    {
      dd = -dd;
    }
    *lr  = (-bb + dd) / aa;
    *ssr = 0.0;
    if(*lr != 0.0)
    {
      *ssr = (cc / *lr) / aa;
    }
    *ssi = 0.0;
    *li  = 0.0;
  }
}
