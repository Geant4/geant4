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

/*
 * Implements calculation of the Fermi density effect as per the method
 * described in:
 *
 *   R. M. Sternheimer, M. J. Berger, and S. M. Seltzer. Density
 *   effect for the ionization loss of charged particles in various sub-
 *   stances. Atom. Data Nucl. Data Tabl., 30:261, 1984.
 *
 * Which (among other Sternheimer references) builds on:
 *
 *   R. M. Sternheimer. The density effect for ionization loss in
 *   materials. Phys. Rev., 88:851­859, 1952.
 *
 * The returned values of delta are directly from the Sternheimer calculation,
 * and not Sternheimer's popular three-part approximate parameterization
 * introduced in the same paper.
 *
 * Author: Matthew Strait <straitm@umn.edu> 2019
 */

#include "G4ios.hh"
#include "G4DensityEffectCalc.hh"
#include "G4NistManager.hh"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

/* Newton's method for finding roots.  Adapted from G4PolynominalSolver, but
 * without the assumption that the input is a polynomial.  Also, here we
 * always expect the roots to be positive, so return -1 as an error value. */
static double newton(const double start,
                     double(*Function)(double, G4DensityEffectCalcData *),
                     double(*Derivative)(double, G4DensityEffectCalcData *),
                     G4DensityEffectCalcData * par)
{
  const int maxIter = 100;

  int nbad = 0, ngood = 0;

  double Lambda = start;

  while(true){
    const double Value = Function(Lambda, par);
    const double Gradient = Derivative(Lambda, par);
    Lambda -= Value/Gradient;

    if(std::fabs((Value/Gradient)/Lambda) <= 1e-12) {
      ngood++;
      if(ngood == 2) return Lambda;
    } else {
      nbad++;
      if(nbad > maxIter || isnan(Value) || isinf(Value)) return -1;
    }
  }
}

/* Return the derivative of the equation used
 * to solve for the Sternheimer parameter rho. */
static double dfrho(double rho, G4DensityEffectCalcData * p)
{
  double ans = 0;
  for(int i = 0; i < p->nlev-1; i++)
    if(p->sternf[i] > 0)
      ans += p->sternf[i] * pow(p->levE[i], 2) * rho /
        (pow(p->levE[i] * rho, 2) + 2./3. * p->sternf[i] * pow(p->plasmaE, 2));

  return ans;
}

/* Return the functional value for the equation used
 * to solve for the Sternheimer parameter rho. */
static double frho(double rho, G4DensityEffectCalcData * p)
{
  double ans = 0;
  for(int i = 0; i < p->nlev-1; i++)
    if(p->sternf[i] > 0)
      ans += p->sternf[i] * log(pow(p->levE[i]*rho, 2) +
        2./3. * p->sternf[i]*pow(p->plasmaE, 2));

  ans *= 0.5; // pulled out of loop for efficiency

  if(p->sternf[p->nlev-1] > 0)
    ans += p->sternf[p->nlev-1] * log(p->plasmaE * sqrt(p->sternf[p->nlev-1]));

  ans -= log(p->meanexcite);

  return ans;
}

/* Return the derivative for the equation used to
 * solve for the Sternheimer parameter l, called 'L' here. */
static double dell(double L, G4DensityEffectCalcData * p)
{
  double ans = 0;
  for(int i = 0; i < p->nlev; i++)
    if(p->sternf[i] > 0 && (p->sternEbar[i] > 0 || L != 0))
      ans += p->sternf[i]/pow(pow(p->sternEbar[i], 2) + L*L, 2);

  ans *= -2*L; // pulled out of the loop for efficiency

  return ans;
}

/* Return the functional value for the equation used to
 * solve for the Sternheimer parameter l, called 'L' here. */
static double ell(double L, G4DensityEffectCalcData * p)
{
  double ans = 0;
  for(int i = 0; i < p->nlev; i++)
    if(p->sternf[i] > 0 && (p->sternEbar[i] > 0 || L != 0))
      ans += p->sternf[i]/(pow(p->sternEbar[i], 2) + L*L);

  ans -= pow(10, - 2 * p->sternx);

  return ans;
}

/**
 * Given the Sternheimer parameter l^2 (called 'sternL' here), and that
 * the l_i and adjusted energies have been found with SetupFermiDeltaCalc(),
 * return the value of delta.  Helper function for DoFermiDeltaCalc().
 */
static double delta_once_solved(G4DensityEffectCalcData * p,
                                const double sternL)
{
  double ans = 0;
  for(int i = 0; i < p->nlev; i++)
    if(p->sternf[i] > 0)
      ans += p->sternf[i] *
        log((pow(p->sternl[i], 2) + pow(sternL, 2))/pow(p->sternl[i], 2));

  ans -= pow(sternL, 2)/(1 + pow(10, 2 * p->sternx));
  return ans;
}

/**
 * Calculate the Sternheimer adjusted energy levels and parameters l_i given
 * the Sternheimer parameter rho.  Helper function for SetupFermiDeltaCalc().
 */
static void calc_Ebarl(G4DensityEffectCalcData * par,
                       const double sternrho)
{
  par->sternl    = (double *)malloc(sizeof(double)*par->nlev);
  par->sternEbar = (double *)malloc(sizeof(double)*par->nlev);
  for(int i = 0; i < par->nlev; i++){
    par->sternEbar[i] = par->levE[i] * sternrho/par->plasmaE;
    par->sternl[i] = i < par->nlev-1?
      sqrt(pow(par->sternEbar[i], 2) + 2./3. * par->sternf[i])
      : sqrt(par->sternf[i]);
  }
}

/**
 * If the "oscillator strengths" par->sternf do not add up to 1,
 * normalize them so that they do.
 */
static void normalize_sternf(G4DensityEffectCalcData * par)
{
  double sum = 0;
  for(int i = 0; i < par->nlev; i++) sum += par->sternf[i];
  for(int i = 0; i < par->nlev; i++) par->sternf[i] /= sum;
}

/**
 * At sufficiently high energy, the density effect takes on a simple
 * form.  Also at sufficiently high energy, we run into numerical problems
 * doing the full calculation.  In that case, return this.
 */
static double delta_limiting_case(G4DensityEffectCalcData * par,
                                  const double sternx)
{
  return 2 * (log(par->plasmaE/par->meanexcite) + log(10.)*sternx) - 1;
}


/********************/
/* Public functions */
/********************/

bool SetupFermiDeltaCalc(G4DensityEffectCalcData * par)
{
  normalize_sternf(par);

  const double sternrho = newton(1.5, frho, dfrho, par);

  // Negative values, and values much larger than unity are non-physical.
  // Values between zero and one are also suspect, but not as clearly wrong.
  if(sternrho <= 0 || sternrho > 100){
    if(G4NistManager::Instance()->GetVerbose() > -1){
      G4cerr <<
      "*********************************************************************\n"
      "* Warning: Could not solve for Sternheimer rho. Probably you have a *\n"
      "* mean ionization energy which is incompatible with your            *\n"
      "* distribution of energy levels, or an unusually dense material.    *\n"
      "* Falling back to parameterization.                                 *\n"
      "*********************************************************************"
      << G4endl
      << "Number of levels: " << par->nlev << G4endl
      << "Mean ionization energy: " << par->meanexcite << "eV" << G4endl
      << "Plasma energy: " << par->plasmaE << "eV" << G4endl;
      for(int i = 0; i < par->nlev; i++)
        G4cerr << "Level " << i
               << ": strength " << par->sternf[i]
               << ": energy " << par->levE[i] << "eV" << G4endl;
    }
    return false;
  }

  calc_Ebarl(par, sternrho);
  return true;
}

double DoFermiDeltaCalc(G4DensityEffectCalcData * par,
                        const double sternx)
{
  // Above beta*gamma of 10^10, the exact treatment is within machine
  // precision of the limiting case, for ordinary solids, at least. The
  // convergence goes up as the density goes down, but even in a pretty
  // hard vacuum it converges by 10^20. Also, it's hard to imagine how
  // this energy is relevant (x = 20 -> 10^19 GeV for muons). So this
  // is mostly not here for physical reasons, but rather to avoid ugly
  // discontinuities in the return value.
  if(sternx > 20) return delta_limiting_case(par, sternx);

  par->sternx = sternx;

  // The derivative of the function we are solving for is strictly
  // negative for positive (physical) values, so if the value at
  // zero is less than zero, it has no solution, and there is no
  // density effect in the Sternheimer "exact" treatment (which is
  // still an approximation).
  if(ell(0, par) <= 0) return 0;

  // Increase initial guess until it converges.
  int startLi;
  for(startLi = -10; startLi < 30; startLi++){
    const double sternL = newton(pow(2, startLi), &ell, &dell, par);
    if(sternL != -1)
      return delta_once_solved(par, sternL);
  }

  return -1; // Signal the caller to use the Sternheimer approximation,
             // because we have been unable to solve the exact form.
}
