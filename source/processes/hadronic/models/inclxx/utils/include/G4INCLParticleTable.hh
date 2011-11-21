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

#ifndef G4INCLParticleTable_hh
#define G4INCLParticleTable_hh 1

#include <string>

#include "G4INCLParticleType.hh"
#include "G4INCLIParticleDataSource.hh"
#include "G4INCLLogger.hh"

namespace G4INCL {
  class ParticleTable {
  public:
    /// Initialize the particle table and use the default particle
    /// data source (PDS)
    static void initialize();

    static void initialize(IParticleDataSource*);

    /** Delete the particle table */
    static void deletePDS() {
      delete pds;
    }

    /// Get the isospin of a particle
    static G4int getIsospin(const ParticleType t);

    /// Get the native INCL name of the particle
    static std::string getName(const ParticleType t);

    /// Get the native INCL name of the ion
    static std::string getName(const G4int A, const G4int Z);

    /// \brief Convert a string to particle type.
    static ParticleType getParticleType(const std::string &pS);

    // Get the name of the particle from the particle data source
    // implementation (e.g. the Geant4 particle name)
    static std::string getPDSName(const ParticleType t);
    static std::string getPDSName(const G4int A, const G4int Z);

    // Get particle mass (in MeV)
    static G4double getMass(const G4INCL::ParticleType t);
    static G4double getMass(const G4int A, const G4int Z);

    /// \brief Get mass number from particle type
    static G4int getMassNumber(const ParticleType t) {
      switch(t) {
        case Proton:
        case Neutron:
        case DeltaPlusPlus:
        case DeltaPlus:
        case DeltaZero:
        case DeltaMinus:
          return 1;
          break;
        case PiPlus:
        case PiMinus:
        case PiZero:
          return 0;
          break;
        default:
          FATAL("Can't determine mass number for particle type " << t << std::endl);
          std::abort();
          return 0;
          break;
      }
    }

    /// \brief Get charge number from particle type
    static G4int getChargeNumber(const ParticleType t) {
      switch(t) {
        case DeltaPlusPlus:
          return 2;
          break;
        case Proton:
        case DeltaPlus:
        case PiPlus:
          return 1;
          break;
        case Neutron:
        case DeltaZero:
        case PiZero:
          return 0;
          break;
        case DeltaMinus:
        case PiMinus:
          return -1;
          break;
        default:
          FATAL("Can't determine charge number for particle type " << t << std::endl);
          std::abort();
          return 0;
          break;
      }
    }

    /// \brief Get RMS radius for clusters
    static G4double getClusterRMS(const G4int A, const G4int Z);

    static G4double getNuclearRadius(const G4int A, const G4int Z);
    static G4double getMaximumNuclearRadius(const G4int A, const G4int Z);
    static G4double getSurfaceDiffuseness(const G4int A, const G4int Z);

    /// \brief Get the separation energy for a particle type.
    static G4double getSeparationEnergy(const ParticleType t) {
      if(t==Proton || t==DeltaPlusPlus || t==DeltaPlus)
        return protonSeparationEnergy;
      else if(t==Neutron || t==DeltaZero || t==DeltaMinus)
        return neutronSeparationEnergy;
      else {
        ERROR("ParticleTable::getSeparationEnergy : Unknown particle type." << std::endl);
        return 0.0;
      }
    }

    const static G4double effectiveNucleonMass;
    const static G4double effectiveNucleonMass2;
    const static G4double effectiveDeltaMass;
    const static G4double effectivePionMass;
    const static G4double effectiveDeltaDecayThreshold;

    static const G4int maxClusterMass = 12;
    static const G4int maxClusterCharge = 8;

    const static G4int clusterTableZSize = ParticleTable::maxClusterCharge+1;
    const static G4int clusterTableASize = ParticleTable::maxClusterMass+1;
    const static G4double binding[clusterTableZSize][clusterTableASize];
    const static G4double rmsc[clusterTableZSize][clusterTableASize];
    const static G4double clusterPosFact[maxClusterMass+1];
    const static G4double clusterPosFact2[maxClusterMass+1];
    const static G4int clusterZMin[maxClusterMass+1]; // Lower limit of Z for cluster of mass A
    const static G4int clusterZMax[maxClusterMass+1]; // Upper limit of Z for cluster of mass A
    const static G4double clusterPhaseSpaceCut[maxClusterMass+1];

    // Enumerator for cluster-decay channels
    enum ClusterDecayType {
      StableCluster,
      NeutronDecay,
      ProtonDecay,
      AlphaDecay,
      TwoProtonDecay,
      TwoNeutronDecay
    };
    const static ClusterDecayType clusterDecayMode[clusterTableZSize][clusterTableASize];

  protected:
    ParticleTable() {};
    ~ParticleTable() {};

  private:
    static IParticleDataSource *pds;
    static G4double protonMass;
    static G4double neutronMass;
    static G4double piPlusMass, piMinusMass, piZeroMass;

    const static G4int mediumNucleiTableSize = 30;
    const static G4double mediumDiffuseness[mediumNucleiTableSize];
    const static G4double pf[mediumNucleiTableSize];
    const static G4double mediumRadius[mediumNucleiTableSize];
    const static G4double rms[mediumNucleiTableSize];
    
    const static G4double neutronSeparationEnergy, protonSeparationEnergy;

  };
}

#endif
