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

#include "G4DensityEffectCalculator.hh"
#include "G4AtomicShells.hh"
#include "G4NistManager.hh"
#include "G4Pow.hh"

static G4Pow * gpow = G4Pow::GetInstance();

const G4int maxWarnings = 20;

G4DensityEffectCalculator::G4DensityEffectCalculator(const G4Material* mat, G4int n)
  : fMaterial(mat), fVerbose(0), fWarnings(0), nlev(n)
{
  fVerbose = std::max(fVerbose, G4NistManager::Instance()->GetVerbose());

  sternf    = new G4double [nlev];
  levE      = new G4double [nlev];
  sternl    = new G4double [nlev];
  sternEbar = new G4double [nlev];
  for(G4int i=0; i<nlev; ++i) {
    sternf[i]    = 0.0;
    levE[i]      = 0.0;
    sternl[i]    = 0.0;
    sternEbar[i] = 0.0;
  }

  fConductivity = sternx = 0.0;
  G4bool conductor = (fMaterial->GetFreeElectronDensity() > 0.0);

  G4int sh = 0;
  G4double sum = 0.;
  const G4double tot = fMaterial->GetTotNbOfAtomsPerVolume();
  for(size_t j = 0; j < fMaterial->GetNumberOfElements(); ++j) {
    // The last subshell is considered to contain the conduction
    // electrons. Sternheimer 1984 says "the lowest chemical valance of
    // the element" is used to set the number of conduction electrons.
    // I'm not sure if that means the highest subshell or the whole
    // shell, but in any case, he also says that the choice is arbitrary
    // and offers a possible alternative. This is one of the sources of
    // uncertainty in the model.
    const G4double frac = fMaterial->GetVecNbOfAtomsPerVolume()[j]/tot;
    const G4int Z = fMaterial->GetElement((G4int)j)->GetZasInt();
    const G4int nshell = G4AtomicShells::GetNumberOfShells(Z);
    for(G4int i = 0; i < nshell; ++i) {
      // For conductors, put *all* top shell electrons into the conduction
      // band, regardless of element.
      const G4double xx = frac*G4AtomicShells::GetNumberOfElectrons(Z, i);
      if(i < nshell-1 || !conductor) {
	sternf[sh] += xx;
      } else {
        fConductivity += xx;
      }
      levE[sh] = G4AtomicShells::GetBindingEnergy(Z, i)/CLHEP::eV;
      ++sh;
    }
  }
  for(G4int i=0; i<nlev; ++i) {
    sum += sternf[i];
  }
  sum += fConductivity;

  const G4double invsum = (sum > 0.0) ? 1./sum : 0.0;
  for(G4int i=0; i<nlev; ++i) {
    sternf[i] *= invsum;
  }
  fConductivity *= invsum;  
  plasmaE = fMaterial->GetIonisation()->GetPlasmaEnergy()/CLHEP::eV;
  meanexcite = fMaterial->GetIonisation()->GetMeanExcitationEnergy()/CLHEP::eV;
}

G4DensityEffectCalculator::~G4DensityEffectCalculator()
{
  delete [] sternf;
  delete [] levE;
  delete [] sternl;
  delete [] sternEbar;
}

G4double G4DensityEffectCalculator::ComputeDensityCorrection(G4double x)
{
  if(fVerbose > 1) {
    G4cout << "G4DensityEffectCalculator::ComputeDensityCorrection for " 
           << fMaterial->GetName() << ", x= " << x << G4endl;
  }
  const G4double approx = fMaterial->GetIonisation()->GetDensityCorrection(x);
  const G4double exact  = FermiDeltaCalculation(x);

  if(fVerbose > 1) {
    G4cout << "   Delta: computed= " << exact 
	   << ", parametrized= " << approx << G4endl;
  }
  if(approx >= 0. && exact < 0.) {
    if(fVerbose > 0) {
      ++fWarnings;
      if(fWarnings < maxWarnings) {
        G4ExceptionDescription ed; 
        ed << "Sternheimer fit failed for " << fMaterial->GetName() 
	   << ", x = " << x << ": Delta exact= "
	   << exact << ", approx= " << approx;
	G4Exception("G4DensityEffectCalculator::DensityCorrection", "mat008", 
                    JustWarning, ed);
      }
    }
    return approx;
  }
  // Fall back to approx if exact and approx are very different, under the
  // assumption that this means the exact calculation has gone haywire
  // somehow, with the exception of the case where approx is negative.  I
  // have seen this clearly-wrong result occur for substances with extremely
  // low density (1e-25 g/cc).
  if(approx >= 0. && std::abs(exact - approx) > 1.) {
    if(fVerbose > 0) {
      ++fWarnings;
      if(fWarnings < maxWarnings) {
        G4ExceptionDescription ed; 
        ed << "Sternheimer exact= " << exact << " and approx= "
           << approx << " are too different for "
           << fMaterial->GetName() << ", x = " << x;
	G4Exception("G4DensityEffectCalculator::DensityCorrection", "mat008", 
                    JustWarning, ed);
      }
    }
    return approx;
  }
  return exact;
}

G4double G4DensityEffectCalculator::FermiDeltaCalculation(G4double x)
{
  // Above beta*gamma of 10^10, the exact treatment is within machine
  // precision of the limiting case, for ordinary solids, at least. The
  // convergence goes up as the density goes down, but even in a pretty
  // hard vacuum it converges by 10^20. Also, it's hard to imagine how
  // this energy is relevant (x = 20 -> 10^19 GeV for muons). So this
  // is mostly not here for physical reasons, but rather to avoid ugly
  // discontinuities in the return value.
  if(x > 20.) { return -1.; }

  sternx = x;
  G4double sternrho = Newton(1.5, true);

  // Negative values, and values much larger than unity are non-physical.
  // Values between zero and one are also suspect, but not as clearly wrong.
  if(sternrho <= 0. || sternrho > 100.) {
    if(fVerbose > 0) {
      ++fWarnings;
      if(fWarnings < maxWarnings) {
        G4ExceptionDescription ed; 
        ed << "Sternheimer computation failed for " << fMaterial->GetName() 
	   << ", x = " << x << ":\n"
	   << "Could not solve for Sternheimer rho. Probably you have a \n"
	   << "mean ionization energy which is incompatible with your\n"
	   << "distribution of energy levels, or an unusually dense material.\n"
	   << "Number of levels: " << nlev 
	   << " Mean ionization energy(eV): " << meanexcite 
	   << " Plasma energy(eV): " << plasmaE << "\n";
	for(G4int i = 0; i < nlev; ++i) {
          ed << "Level " << i << ": strength " << sternf[i]
             << ": energy(eV)= " << levE[i] << "\n";
	}
	G4Exception("G4DensityEffectCalculator::SetupFermiDeltaCalc", "mat008", 
                    JustWarning, ed);
      }
    }
    return -1.;
  }

  // Calculate the Sternheimer adjusted energy levels and parameters l_i given
  // the Sternheimer parameter rho.
  for(G4int i=0; i<nlev; ++i) {
    sternEbar[i] = levE[i] * (sternrho/plasmaE);
    sternl[i] = std::sqrt(gpow->powN(sternEbar[i], 2) + (2./3.)*sternf[i]);
  }
  // The derivative of the function we are solving for is strictly
  // negative for positive (physical) values, so if the value at
  // zero is less than zero, it has no solution, and there is no
  // density effect in the Sternheimer "exact" treatment (which is
  // still an approximation).
  //
  // For conductors, this test is not needed, because Ell(L) contains
  // the term fConductivity/(L*L), so the value at L=0 is always
  // positive infinity. In the code we don't return inf, though, but
  // rather set that term to zero, which means that if this test were
  // used, it would give the wrong result for some materials.
  if(fConductivity == 0 && Ell(0) <= 0)
  {
    return 0;
  }

  // Attempt to find the root from 40 starting points evenly distributed
  // in log space.  Trying a single starting point is not sufficient for
  // convergence in most cases.
  for(G4int startLi = -10; startLi < 30; ++startLi){
    const G4double sternL = Newton(gpow->powN(2, startLi), false);
    if(sternL != -1.) {
      return DeltaOnceSolved(sternL);
    }
  }
  return -1.; // Signal the caller to use the Sternheimer approximation,
              // because we have been unable to solve the exact form.
}

/* Newton's method for finding roots.  Adapted from G4PolynominalSolver, but
 * without the assumption that the input is a polynomial.  Also, here we
 * always expect the roots to be positive, so return -1 as an error value. */
G4double G4DensityEffectCalculator::Newton(G4double start, G4bool first)
{
  const G4int maxIter = 100;
  G4int nbad = 0, ngood = 0;

  G4double lambda(start), value(0.), dvalue(0.);

  if(fVerbose > 2) {
    G4cout << "G4DensityEffectCalculator::Newton: strat= " << start 
	   << " type: " << first << G4endl; 
  }
  while(true) {
    if(first) {
      value = FRho(lambda);
      dvalue = DFRho(lambda);
    } else {
      value = Ell(lambda);
      dvalue = DEll(lambda);
    }
    if(dvalue == 0.0) { break; }
    const G4double del = value/dvalue;
    lambda -= del;

    const G4double eps = std::abs(del/lambda);
    if(eps <= 1.e-12) {
      ++ngood;
      if(ngood == 2) { 
	if(fVerbose > 2) {
	  G4cout << "  Converged with result= " << lambda << G4endl; 
	}
	return lambda; 
      }
    } else {
      ++nbad;
    }
    if(nbad > maxIter || std::isnan(value) || std::isinf(value)) { break; }
  }
  if(fVerbose > 2) {
    G4cout << "  Failed to converge last value= " << value 
	   << " dvalue= " << dvalue << " lambda= " << lambda << G4endl; 
  }
  return -1.;
}

/* Return the derivative of the equation used
 * to solve for the Sternheimer parameter rho. */
G4double G4DensityEffectCalculator::DFRho(G4double rho)
{
  G4double ans = 0.0;
  for(G4int i = 0; i < nlev; ++i) {
    if(sternf[i] > 0.) {
      ans += sternf[i] * gpow->powN(levE[i], 2) * rho /
        (gpow->powN(levE[i] * rho, 2)
         + 2./3. * sternf[i] * gpow->powN(plasmaE, 2));
    }
  }
  return ans;
}

/* Return the functional value for the equation used
 * to solve for the Sternheimer parameter rho. */
G4double G4DensityEffectCalculator::FRho(G4double rho)
{
  G4double ans = 0.0;
  for(G4int i = 0; i<nlev; ++i) {
    if(sternf[i] > 0.) {
      ans += sternf[i] * G4Log(gpow->powN(levE[i]*rho, 2) +
        2./3. * sternf[i]*gpow->powN(plasmaE, 2));
    }    
  }
  ans *= 0.5; // pulled out of loop for efficiency

  if(fConductivity > 0.) {
    ans += fConductivity * G4Log(plasmaE * std::sqrt(fConductivity));
  }
  ans -= G4Log(meanexcite);
  return ans;
}

/* Return the derivative for the equation used to
 * solve for the Sternheimer parameter l, called 'L' here. */
G4double G4DensityEffectCalculator::DEll(G4double L)
{
  G4double ans = 0.;
  for(G4int i=0; i<nlev; ++i) {
    if(sternf[i] > 0 && (sternEbar[i] > 0. || L != 0.)) {
      const G4double y = gpow->powN(sternEbar[i], 2);
      ans += sternf[i]/gpow->powN(y + L*L, 2);
    }
  }
  ans += fConductivity/gpow->powN(L*L, 2);
  ans *= (-2*L); // pulled out of the loop for efficiency
  return ans;
}

/* Return the functional value for the equation used to
 * solve for the Sternheimer parameter l, called 'L' here. */
G4double G4DensityEffectCalculator::Ell(G4double L)
{
  G4double ans = 0.;
  for(G4int i=0; i<nlev; ++i) {
    if(sternf[i] > 0. && (sternEbar[i] > 0. || L != 0.)) {
      ans += sternf[i]/(gpow->powN(sternEbar[i], 2) + L*L);
    }
  }
  if(fConductivity > 0. && L != 0.) {
    ans += fConductivity/(L*L);
  }  
  ans -= gpow->powZ(10, -2 * sternx);
  return ans;
}

/**
 * Given the Sternheimer parameter l (called 'sternL' here), and that
 * the l_i and adjusted energies have been found with SetupFermiDeltaCalc(),
 * return the value of delta.  Helper function for DoFermiDeltaCalc().
 */
G4double G4DensityEffectCalculator::DeltaOnceSolved(G4double sternL)
{
  G4double ans = 0.;
  for(G4int i=0; i<nlev; ++i) {
    if(sternf[i] > 0.) {
      ans += sternf[i] * G4Log((gpow->powN(sternl[i], 2)
              + gpow->powN(sternL, 2))/gpow->powN(sternl[i], 2));
    }
  }
  // sternl for the conduction electrons is sqrt(fConductivity), with
  // no factor of 2./3 as with the other levels.
  if(fConductivity > 0) {
    ans += fConductivity * G4Log((fConductivity
                  + gpow->powN(sternL, 2))/fConductivity);
  }
  ans -= gpow->powN(sternL, 2)/(1 + gpow->powZ(10, 2 * sternx));
  return ans;
}
