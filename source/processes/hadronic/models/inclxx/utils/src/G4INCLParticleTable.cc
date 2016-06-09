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
// Pekka Kaitaniemi, CEA and Helsinki Institute of Physics
// Davide Mancusi, CEA
// Alain Boudard, CEA
// Sylvie Leray, CEA
// Joseph Cugnon, University of Liege
//
// INCL++ revision: v5.0_rc3
//
#define INCLXX_IN_GEANT4_MODE 1

#include "globals.hh"

#include "G4INCLParticleTable.hh"
#include "G4INCLStandaloneParticleDataSource.hh"
#include <algorithm>
// Assert is banned for now, for G4 compatibility reasons...
//#include <cassert>

namespace G4INCL {

  IParticleDataSource* ParticleTable::pds = 0;
  G4double ParticleTable::protonMass = 0.0;
  G4double ParticleTable::neutronMass = 0.0;
  G4double ParticleTable::piPlusMass = 0.0;
  G4double ParticleTable::piMinusMass = 0.0;
  G4double ParticleTable::piZeroMass = 0.0;
  const G4double ParticleTable::mediumDiffuseness[mediumNucleiTableSize] =
  {0.0,0.0,0.0,0.0,0.0,1.78,1.77,1.77,1.77,1.71,
    1.69,1.69,1.635,1.730,1.81,1.833,1.798,
    1.841,0.567,0.571, 0.560,0.549,0.550,0.551,
    0.580,0.575,0.569,0.537,0.0,0.0};
  const G4double ParticleTable::pf[mediumNucleiTableSize] = {0.0, 0.0, 0.0, 0.0, 0.0, 77.0, 110.0, 110.0, 153.0};
  const G4double ParticleTable::mediumRadius[mediumNucleiTableSize] =
  {0.0,0.0,0.0,0.0,0.0,0.334,0.327,0.479,0.631,0.838,
    0.811,1.07,1.403,1.335,1.25,1.544,1.498,1.513,
    2.58,2.77, 2.775,2.78,2.88,2.98,3.22,3.03,2.84,
    3.14,0.0,0.0};
  const G4double ParticleTable::rms[mediumNucleiTableSize] = {0.0, 0.0, 0.0, 0.0, 0.0, 2.10, 1.80, 1.80, 1.63};

  const G4double ParticleTable::binding[clusterTableZSize][clusterTableASize] =
  {
    {0.00, 0.00, 0.00, 0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00},
    {0.00, 0.00, 2.22, 8.48,  5.60,  6.68,  5.76,  6.58,  0.00,  0.00,  0.00,  0.00,  0.00},
    {0.00, 0.00, 0.00, 7.72, 28.30, 27.41, 29.27, 28.83, 31.40, 31.29, 30.34,  0.00,  0.00},
    {0.00, 0.00, 0.00, 0.00,  4.60, 26.33, 31.99, 39.25, 41.28, 45.34, 45.31, 45.71, 44.40},
    {0.00, 0.00, 0.00, 0.00,  0.00,  0.00, 26.92, 37.60, 56.50, 58.16, 64.98, 65.48, 68.65},
    {0.00, 0.00, 0.00, 0.00,  0.00,  0.00,  0.90, 24.71, 37.74, 56.31, 64.75, 76.20, 79.58},
    {0.00, 0.00, 0.00, 0.00,  0.00,  0.00,  0.00,  0.00, 24.78, 39.04, 60.32, 73.44, 92.16},
    {0.00, 0.00, 0.00, 0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00, 36.40, 59.00, 74.04},
    {0.00, 0.00, 0.00, 0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00, 58.549}
  };

  const G4double ParticleTable::rmsc[clusterTableZSize][clusterTableASize] =
  {
    {0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00},
    {0.00, 0.00, 2.10, 1.80, 1.70, 1.83, 2.60, 2.50, 0.00, 0.00, 0.00, 0.00, 0.00},
    {0.00, 0.00, 0.00, 1.80, 1.68, 1.70, 2.60, 2.50, 2.50, 2.50, 2.50, 0.00, 0.00},
    {0.00, 0.00, 0.00, 0.00, 1.70, 1.83, 2.56, 2.40, 2.50, 2.50, 2.50, 2.50, 2.50},
    {0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 2.60, 2.50, 2.50, 2.51, 2.50, 2.50, 2.50},
    {0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 2.50, 2.50, 2.50, 2.50, 2.45, 2.40, 2.50},
    {0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 2.50, 2.50, 2.50, 2.50, 2.47},
    {0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 2.50, 2.50, 2.50},
    {0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 2.50}
  };

  /* Table for cluster decays
   * Definition of "Stable": halflife > 1 ms
   *
   * These table includes decay data for clusters that INCL presently does
   * not produce. It can't hurt.
   *
   * Unbound nuclides (like C-6) are marked as stable, but should never be
   * produced by INCL. If you find them in the output, something is fishy.
   */
  const ParticleTable::ClusterDecayType ParticleTable::clusterDecayMode[clusterTableZSize][clusterTableASize] =
  {
    {StableCluster, StableCluster, StableCluster, StableCluster, StableCluster, StableCluster, StableCluster, StableCluster, StableCluster, StableCluster, StableCluster, StableCluster, StableCluster},
    {StableCluster, StableCluster, StableCluster, StableCluster, NeutronDecay, TwoNeutronDecay, NeutronDecay, TwoNeutronDecay, StableCluster, StableCluster, StableCluster, StableCluster, StableCluster},
    {StableCluster, StableCluster, StableCluster, StableCluster, StableCluster, NeutronDecay, StableCluster, NeutronDecay, StableCluster, NeutronDecay, NeutronDecay, StableCluster, StableCluster},
    {StableCluster, StableCluster, StableCluster, StableCluster, ProtonDecay, ProtonDecay, StableCluster, StableCluster, StableCluster, StableCluster, NeutronDecay, StableCluster, NeutronDecay},
    {StableCluster, StableCluster, StableCluster, StableCluster, StableCluster, StableCluster, TwoProtonDecay, StableCluster, AlphaDecay, StableCluster, StableCluster, StableCluster, StableCluster},
    {StableCluster, StableCluster, StableCluster, StableCluster, StableCluster, StableCluster, TwoProtonDecay, ProtonDecay, StableCluster, ProtonDecay, StableCluster, StableCluster, StableCluster},
    {StableCluster, StableCluster, StableCluster, StableCluster, StableCluster, StableCluster, StableCluster, StableCluster, TwoProtonDecay, StableCluster, StableCluster, StableCluster, StableCluster},
    {StableCluster, StableCluster, StableCluster, StableCluster, StableCluster, StableCluster, StableCluster, StableCluster, StableCluster, StableCluster, ProtonDecay, ProtonDecay, StableCluster},
    {StableCluster, StableCluster, StableCluster, StableCluster, StableCluster, StableCluster, StableCluster, StableCluster, StableCluster, StableCluster, StableCluster, StableCluster, ProtonDecay}
  };

  /**
   * Precomputed factor 1.0/A
   */
  const G4double ParticleTable::clusterPosFact[maxClusterMass+1] = {0.0, 1.0, 0.5,
						      0.33333, 0.25,
						      0.2, 0.16667,
						      0.14286, 0.125,
						      0.11111, 0.1,
						      0.09091, 0.083333};

  /**
   * Precomputed factor (1.0/A)^2
   */
  const G4double ParticleTable::clusterPosFact2[maxClusterMass+1] = {0.0, 1.0, 0.25,
						       0.11111, 0.0625,
						       0.04, 0.0277778,
						       0.020408, 0.015625,
						       0.012346, 0.01,
						       0.0082645, 0.0069444};

  const G4int ParticleTable::clusterZMin[maxClusterMass+1] = {0, 0, 1, 1, 1, 1, 1, 1, 2, 2, 2, 3, 3};
  const G4int ParticleTable::clusterZMax[maxClusterMass+1] = {0, 0, 1, 2, 3, 3, 5, 5, 6, 6, 7, 7, 8};
  const G4double ParticleTable::clusterPhaseSpaceCut[maxClusterMass+1] = {0.0, 70000.0, 180000.0,
							    90000.0, 90000.0,
							    128941.0 ,145607.0,
							    161365.0, 176389.0,
							    190798.0, 204681.0,
							    218109.0, 231135.0};
  const G4double ParticleTable::effectiveNucleonMass = 938.2796;
  const G4double ParticleTable::effectiveNucleonMass2 = 8.8036860777616e5;
  const G4double ParticleTable::effectiveDeltaMass = 1232.0;
  const G4double ParticleTable::effectivePionMass = 138.0;
  const G4double ParticleTable::effectiveDeltaDecayThreshold =
    ParticleTable::effectiveNucleonMass + ParticleTable::effectivePionMass
    + 2.0;
  const G4double ParticleTable::protonSeparationEnergy = 6.83;
  const G4double ParticleTable::neutronSeparationEnergy = ParticleTable::protonSeparationEnergy;

  void ParticleTable::initialize() {
    ParticleTable::pds = new StandaloneParticleDataSource;
    protonMass = pds->getMass(Proton);
    neutronMass = pds->getMass(Neutron);
    piPlusMass = pds->getMass(PiPlus);
    piMinusMass = pds->getMass(PiMinus);
    piZeroMass = pds->getMass(PiZero);
  }

  void ParticleTable::initialize(IParticleDataSource *p) {
    ParticleTable::pds = p;
  }

  G4int ParticleTable::getIsospin(const ParticleType t) {
    // Actually this is the 3rd component of isospin (I_z) multiplied by 2!
    if(t == Proton) {
      return 1;
    } else if(t == Neutron) {
      return -1;
    } else if(t == PiPlus) {
      return 2;
    } else if(t == PiMinus) {
      return -2;
    } else if(t == PiZero) {
      return 0;
    } else if(t == DeltaPlusPlus) {
      return 3;
    } else if(t == DeltaPlus) {
      return 1;
    } else if(t == DeltaZero) {
      return -1;
    } else if(t == DeltaMinus) {
      return -3;
    }

    ERROR("Requested isospin of an unknown particle!");
    return -10; // Unknown
  }

  std::string ParticleTable::getName(const ParticleType p) {
    if(p == G4INCL::Proton) {
      return std::string("proton");
    } else if(p == G4INCL::Neutron) {
      return std::string("neutron");
    } else if(p == G4INCL::DeltaPlusPlus) {
      return std::string("delta++");
    } else if(p == G4INCL::DeltaPlus) {
      return std::string("delta+");
    } else if(p == G4INCL::DeltaZero) {
      return std::string("delta0");
    } else if(p == G4INCL::DeltaMinus) {
      return std::string("delta-");
    } else if(p == G4INCL::PiPlus) {
      return std::string("pi+");
    } else if(p == G4INCL::PiZero) {
      return std::string("pi0");
    } else if(p == G4INCL::PiMinus) {
      return std::string("pi-");
    }
    return std::string("unknown");
  }

  ParticleType ParticleTable::getParticleType(const std::string &pS) {
    // Normalise the string to lower case
    std::string pSNorm = pS;
    std::transform(pSNorm.begin(), pSNorm.end(), pSNorm.begin(), ::tolower);
    if(pSNorm=="p" || pSNorm=="proton")
      return G4INCL::Proton;
    else if(pSNorm=="n" || pSNorm=="neutron")
      return G4INCL::Neutron;
    else if(pSNorm=="delta++" || pSNorm=="deltaplusplus")
      return G4INCL::DeltaPlusPlus;
    else if(pSNorm=="delta+" || pSNorm=="deltaplus")
      return G4INCL::DeltaPlus;
    else if(pSNorm=="delta0" || pSNorm=="deltazero")
      return G4INCL::DeltaZero;
    else if(pSNorm=="delta-" || pSNorm=="deltaminus")
      return G4INCL::DeltaMinus;
    else if(pSNorm=="pi+" || pSNorm=="pion+" || pSNorm=="piplus" || pSNorm=="pionplus")
      return G4INCL::PiPlus;
    else if(pSNorm=="pi0" || pSNorm=="pion0" || pSNorm=="pizero" || pSNorm=="pionzero")
      return G4INCL::PiZero;
    else if(pSNorm=="pi-" || pSNorm=="pion-" || pSNorm=="piminus" || pSNorm=="pionminus")
      return G4INCL::PiMinus;
    else
      return G4INCL::UnknownParticle;
  }

  G4double ParticleTable::getMass(const ParticleType pt) {
    if(pt == Proton) {
      return protonMass;
    } else if(pt == Neutron) {
      return neutronMass;
    } else if(pt == PiPlus) {
      return piPlusMass;
    } else if(pt == PiMinus) {
      return piMinusMass;
    } else if(pt == PiZero) {
      return piZeroMass;
    } else {
      ERROR("ParticleTable::getMass : Unknown particle type." << std::endl);
      return 0.0;
    }
  }

  G4double ParticleTable::getMass(const G4int A, const G4int Z) {
    //    assert(A>0 && Z>=0 && Z<=A);
    if(A<=clusterTableASize && Z<=clusterTableZSize)
      // should return the correct masses for protons and neutrons, too
      return Z*protonMass + (A-Z)*neutronMass - binding[Z][A];
    else
      return Z*(protonMass - protonSeparationEnergy) + (A-Z)*(neutronMass - neutronSeparationEnergy);
  }

  G4double ParticleTable::getClusterRMS(const G4int A, const G4int Z) {
    //    assert(A>0 && Z>=0 && Z<=A && A<=maxClusterMass && Z<=maxClusterCharge);
    return rmsc[Z][A];
  }

  G4double ParticleTable::getNuclearRadius(const G4int A, const G4int Z) {
    if(A >= 28) {
      return (2.745e-4 * A + 1.063) * std::pow(A, 1.0/3.0);
    } else if(A < 28 && A >= 19) {
      return mediumRadius[A-1];
    } else if(A < 19 && A >= 6) {
      return mediumRadius[A-1];
      //      return 1.581*mediumDiffuseness[A-1]*(2.+5.*mediumRadius[A-1])/(2.+3.*mediumRadius[A-1]);
    } else if(A < 6 && A >= 3) {
      if(A == 3 && Z == 1) {
        return rms[6];
      } else if(A == 3 && Z == 2) {
        return rms[7];
      } else if(A == 4) {
        return rms[8];
      } else {
        ERROR("ParticleTable::getNuclearRadius : No radius for nucleus A = " << A << " Z = " << Z << std::endl);
        return 0.0;
      }
    } else if(A == 2) {
      return rms[5];
    } else {
      ERROR("ParticleTable::getNuclearRadius : No radius for nucleus A = " << A << " Z = " << Z << std::endl);
      return 0.0;
    }
  }

  G4double ParticleTable::getMaximumNuclearRadius(const G4int A, const G4int Z) {
    const G4double XFOISA = 8.0;
    const G4double radius = getNuclearRadius(A,Z);
    const G4double diffuseness = getSurfaceDiffuseness(A,Z);
    if(A >= 28) {
      return radius + XFOISA * diffuseness;
    } else if(A < 28 && A >= 19) {
      return radius + XFOISA * diffuseness;
    } else if(A < 19 && A >= 6) {
      return 5.5 + 0.3 * (G4double(A) - 6.0)/12.0;
    } else if(A >= 2) {
      return ParticleTable::getNuclearRadius(A, Z) + 2.5;
    } else {
      ERROR("ParticleTable::getMaximumNuclearRadius : No maximum radius for nucleus A = " << A << " Z = " << Z << std::endl);
      return 0.0;
    }
  }

  G4double ParticleTable::getSurfaceDiffuseness(const G4int A, const G4int Z) {

    if(A >= 28) {
      return 1.63e-4 * A + 0.510;
    } else if(A < 28 && A >= 19) {
      return mediumDiffuseness[A-1];
    } else if(A < 19 && A >= 6) {
      return mediumDiffuseness[A-1];
    } else if(A < 6 && A >= 2) {
      return 0.57735 * ParticleTable::getNuclearRadius(A, Z);
    } else {
      ERROR("ParticleTable::getSurfaceDiffuseness : No diffuseness for nucleus A = " << A << " Z = " << Z << std::endl);
      return 0.0;
    }
  }
}
