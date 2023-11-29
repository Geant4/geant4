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

/*
 * Interface to calculation of the Fermi density effect as per the method
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

#ifndef G4DensityEffectCalculator_HH
#define G4DensityEffectCalculator_HH

#include "globals.hh"

class G4Material;

class G4DensityEffectCalculator
{
 public:
  G4DensityEffectCalculator(const G4Material*, G4int);
  ~G4DensityEffectCalculator();

  // The Sternheimer 'x' defined as log10(p/m) == log10(beta*gamma).
  G4double ComputeDensityCorrection(G4double x);

 private:
  /*
   * Given a material defined in 'par' with a plasma energy, mean excitation
   * energy, and set of atomic energy levels ("oscillator frequencies") with
   * occupation fractions ("oscillation strengths"), solve for the Sternheimer
   * adjustment factor (Sternheimer 1984 eq 8) and record (into 'par') the values
   * of the adjusted oscillator frequencies and Sternheimer constants l_i.
   * After doing this, 'par' is ready for a calculation of delta for an
   * arbitrary particle energy.  Returns true on success, false on failure.
   */
  G4double FermiDeltaCalculation(G4double x);

  G4double Newton(G4double x0, G4bool first);

  G4double DFRho(G4double);

  G4double FRho(G4double);

  G4double DEll(G4double);

  G4double Ell(G4double);

  G4double DeltaOnceSolved(G4double);

  const G4Material* fMaterial;
  G4int fVerbose{0};
  G4int fWarnings{0};

  // Number of energy levels.  If a single element, this is the number
  // of subshells.  If several elements, this is the sum of the number
  // of subshells.  In principle, could include levels for molecular
  // orbitals or other non-atomic states.  The last level is always
  // the conduction band.  If the material is an insulator, set the
  // oscillator strength for that level to zero and the energy to
  // any value.
  const G4int nlev;

  G4double fConductivity;

  // Current Sternheimer 'x' defined as log10(p/m) == log10(beta*gamma).
  G4double sternx;

  // The plasma energy of the material in eV, which is simply
  // 28.816 sqrt(density Z/A), with density in g/cc.
  G4double plasmaE;

  // The mean excitation energy of the material in eV, i.e. the 'I' in the
  // Bethe energy loss formula.
  G4double meanexcite;

  // Sternheimer's "oscillator strengths", which are simply the fraction
  // of electrons in a given energy level.  For a single element, this is
  // the fraction of electrons in a subshell.  For a compound or mixture,
  // it is weighted by the number fraction of electrons contributed by
  // each element, e.g. for water, oxygen's electrons are given 8/10 of the
  // weight.
  G4double* sternf;

  // Energy levels.  Can be found for free atoms in, e.g., T. A. Carlson.
  // Photoelectron and Auger Spectroscopy. Plenum Press, New York and London,
  // 1985. Available in a convenient form in G4AtomicShells.cc.
  //
  // Sternheimer 1984 implies that the energy level for conduction electrons
  // (the final element of this array) should be set to zero, although the
  // computation could be run with other values.
  G4double* levE;

  /***** Results of intermediate calculations *****/

  // The Sternheimer parameters l_i which appear in Sternheimer 1984 eq(1).
  G4double* sternl;

  // The adjusted energy levels, as found using Sternheimer 1984 eq(8).
  G4double* sternEbar;
};

#endif
