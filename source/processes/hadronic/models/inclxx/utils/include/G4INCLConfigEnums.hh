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
// INCL++ intra-nuclear cascade model
// Alain Boudard, CEA-Saclay, France
// Joseph Cugnon, University of Liege, Belgium
// Jean-Christophe David, CEA-Saclay, France
// Pekka Kaitaniemi, CEA-Saclay, France, and Helsinki Institute of Physics, Finland
// Sylvie Leray, CEA-Saclay, France
// Davide Mancusi, CEA-Saclay, France
//
#define INCLXX_IN_GEANT4_MODE 1

#include "globals.hh"

#ifndef G4INCLConfigEnums_hh
#define G4INCLConfigEnums_hh

namespace G4INCL {

  // Enumerator for Pauli-blocking algorithms
  enum PauliType {
    StatisticalPauli,
    StrictPauli,
    StrictStatisticalPauli,
    GlobalPauli,
    NoPauli
  };

  // Enumerator for Coulomb-distortion algorithms
  enum CoulombType {
    NonRelativisticCoulomb,
    NoCoulomb
  };

  // Enumerator for potential types
  enum PotentialType {
    IsospinEnergySmoothPotential,
    IsospinEnergyPotential,
    IsospinPotential,
    ConstantPotential
  };

  // Enumerator for local-energy types
  enum LocalEnergyType {
    AlwaysLocalEnergy,
    FirstCollisionLocalEnergy,
    NeverLocalEnergy
  };

  // Enumerator for de-excitation types
  enum DeExcitationType {
    DeExcitationNone
#ifdef INCL_DEEXCITATION_ABLAXX
    , DeExcitationABLAv3p
#endif
#ifdef INCL_DEEXCITATION_ABLA07
    , DeExcitationABLA07
#endif
#ifdef INCL_DEEXCITATION_SMM
    , DeExcitationSMM
#endif
#ifdef INCL_DEEXCITATION_GEMINIXX
    , DeExcitationGEMINIXX
#endif
  };

  // Enumerator for cluster-algorithm types
  enum ClusterAlgorithmType {
    IntercomparisonClusterAlgorithm,
    NoClusterAlgorithm
  };

  // Enumerator for separation-energy types
  enum SeparationEnergyType {
    INCLSeparationEnergy,
    RealSeparationEnergy,
    RealForLightSeparationEnergy
  };

  // Enumerator for Fermi-momentum types
  enum FermiMomentumType {
    ConstantFermiMomentum,
    ConstantLightFermiMomentum,
    MassDependentFermiMomentum
  };

  // Enumerator for RNG
  enum RNGType {
    RanecuType,
    Ranecu3Type
  };

  // Enumerator for Cross-Section parametrizations
  enum CrossSectionsType {
    INCL46CrossSections,
    MultiPionsCrossSections,
	TruncatedMultiPionsCrossSections,
	MultiPionsAndResonancesCrossSections,
	StrangenessCrossSections
  };

  // Enumerator for phase-space generator
  enum PhaseSpaceGeneratorType {
    KopylovType,
    RauboldLynchType
  };

  // Enumerator for cascade actions
  enum CascadeActionType {
    DefaultActionType,
    AvatarDumpActionType
  };

}

#endif
