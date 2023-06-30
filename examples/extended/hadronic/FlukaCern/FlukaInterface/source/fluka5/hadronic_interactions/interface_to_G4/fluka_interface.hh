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
// Interface to FLUKA hadron inelastic physics.
//
// XS and final states are accessed independently, 
// to match the G4 needs 
// (G4VCrossSectionDataSet for XS + G4HadronicInteraction for FS).
// 
// FLUKA inelastic hadron-nucleus interactions:
// Hadron-NUCLEON interaction models are based on resonance production and decay below a few GeV, 
// and on the Dual Parton model above. 
// Hadron-NUCLEUS interactions: the PEANUT package includes 
// a detailed Generalised Intra-Nuclear Cascade (GINC) and a preequilibrium stage, 
// followed by equilibrium processes: evaporation, fission, Fermi break-up, gamma deexcitation.
// A. Ferrari and P. Sala, “The Physics of High Energy Reactions,” in Proc. Workshop on Nuclear Reaction Data and Nuclear Reactors Physics, Design and Safety, p. 424, World Scientiﬁc, 1998.
// A. Ferrari and P. Sala, “Nuclear reactions in Monte Carlo codes,” Radiat. Prot. Dosimetry, vol. 99, no. 1-4, pp. 29–38, 2002.
//
//
// NB 1: The user can choose, directly in G4 client code, to:
// - activate coalescence.
// - activate heavy fragments evaporation.
// One can simply pass the relevant booleans 
// as parameters to fluka_interface::initialize(...).
//
// NB 2: This interface also provides support for photonuclear reactions,
// though use in this context has been less extensively debugged.
//
// Author: G.Hugo, 01 August 2022
//
// ***************************************************************************
#ifdef G4_USE_FLUKA
#ifndef FLUKA_INTERFACE_HH
#define FLUKA_INTERFACE_HH


#include <utility>
// G4
#include "globals.hh"


class G4DynamicParticle;
class G4HadProjectile;
class G4Nucleus;
class G4HadFinalState;


namespace fluka_interface {

  void initialize(bool activateCoalescence = false,
                  G4bool activateHeavyFragmentsEvaporation = false);
  G4double computeInelasticScatteringXS(const G4DynamicParticle* projectile,
                                      const G4int targetZ,
                                      const G4int targetA = 0);
  void setNuclearInelasticFinalState(G4HadFinalState* const finalState,
                                     const G4HadProjectile& projectile, 
                                     const G4Nucleus& targetNucleus);

  // HELPERS
  void updateFLUKAProjectileId(G4int& kproj);
  void transformNonSupportedHadrons(G4int& kproj, G4double& ekproj);
  std::pair<G4double, G4double> getKineticEnergyAndMomentum(const G4double ekproj, const G4int kproj);
}


#endif
#endif // G4_USE_FLUKA
