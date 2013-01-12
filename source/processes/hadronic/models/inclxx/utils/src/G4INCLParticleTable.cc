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
#define INCLXX_IN_GEANT4_MODE 1

#include "globals.hh"

#include "G4INCLParticleTable.hh"
#include <algorithm>
// #include <cassert>
#include <cmath>
#include <cctype>
#include <sstream>

#ifdef INCLXX_IN_GEANT4_MODE
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#endif

namespace G4INCL {

  /// \brief Static instance of the NaturalIsotopicAbundances class
  const NaturalIsotopicDistributions *ParticleTable::theNaturalIsotopicDistributions = NULL;

  /// \brief Static pointer to the mass function for nuclei
  ParticleTable::NuclearMassFn ParticleTable::getTableMass;
  /// \brief Static pointer to the mass function for particles
  ParticleTable::ParticleMassFn ParticleTable::getTableParticleMass;
  /// \brief Static pointer to the separation-energy function
  ParticleTable::SeparationEnergyFn ParticleTable::getSeparationEnergy;

  const G4double ParticleTable::theINCLNucleonMass = 938.2796;
  const G4double ParticleTable::theINCLPionMass = 138.0;
  G4double ParticleTable::protonMass = 0.0;
  G4double ParticleTable::neutronMass = 0.0;
  G4double ParticleTable::piPlusMass = 0.0;
  G4double ParticleTable::piMinusMass = 0.0;
  G4double ParticleTable::piZeroMass = 0.0;

  // e^2/(4 pi epsilon_0) [MeV fm]
  const G4double ParticleTable::eSquared = 1.439964;

  // Hard-coded values of the real particle masses (MeV/c^2)
  G4double ParticleTable::theRealProtonMass = 938.27203;
  G4double ParticleTable::theRealNeutronMass = 939.56536;
  G4double ParticleTable::theRealChargedPiMass = 139.57018;
  G4double ParticleTable::theRealPiZeroMass = 134.9766;

  const G4double ParticleTable::mediumDiffuseness[mediumNucleiTableSize] =
  {0.0,0.0,0.0,0.0,0.0,1.78,1.77,1.77,1.77,1.71,
    1.69,1.69,1.635,1.730,1.81,1.833,1.798,
    1.841,0.567,0.571, 0.560,0.549,0.550,0.551,
    0.580,0.575,0.569,0.537,0.0,0.0};
  const G4double ParticleTable::mediumRadius[mediumNucleiTableSize] =
  {0.0,0.0,0.0,0.0,0.0,0.334,0.327,0.479,0.631,0.838,
    0.811,1.07,1.403,1.335,1.25,1.544,1.498,1.513,
    2.58,2.77, 2.775,2.78,2.88,2.98,3.22,3.03,2.84,
    3.14,0.0,0.0};
  const G4double ParticleTable::positionRMS[clusterTableZSize][clusterTableASize] = {
/*     A=   0     1     2     3     4     5     6     7     8     9    10    11    12 */
/* Z=0 */ {0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00},
/* Z=1 */ {0.00, 0.00, 2.10, 1.80, 1.70, 1.83, 2.60, 2.50, 0.00, 0.00, 0.00, 0.00, 0.00},
/* Z=2 */ {0.00, 0.00, 0.00, 1.80, 1.68, 1.70, 2.60, 2.50, 2.50, 2.50, 2.50, 0.00, 0.00},
/* Z=3 */ {0.00, 0.00, 0.00, 0.00, 1.70, 1.83, 2.56, 2.40, 2.50, 2.50, 2.50, 2.50, 2.50},
/* Z=4 */ {0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 2.60, 2.50, 2.50, 2.51, 2.50, 2.50, 2.50},
/* Z=5 */ {0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 2.50, 2.50, 2.50, 2.50, 2.45, 2.40, 2.50},
/* Z=6 */ {0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 2.50, 2.50, 2.50, 2.50, 2.47},
/* Z=7 */ {0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 2.50, 2.50, 2.50},
/* Z=8 */ {0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 2.50}
  };

  const G4double ParticleTable::momentumRMS[clusterTableZSize][clusterTableASize] = {
/*     A=   0     1     2     3     4     5     6     7     8     9    10    11    12 */
/* Z=0 */ {0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00},
/* Z=1 */ {0.00, 0.00, 77.0, 110., 153., 100., 100., 100., 0.00, 0.00, 0.00, 0.00, 0.00},
/* Z=2 */ {0.00, 0.00, 0.00, 110., 153., 100., 100., 100., 100., 100., 100., 0.00, 0.00},
/* Z=3 */ {0.00, 0.00, 0.00, 0.00, 153., 100., 100., 100., 100., 100., 100., 100., 100.},
/* Z=4 */ {0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 100., 100., 100., 100., 100., 100., 100.},
/* Z=5 */ {0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 100., 100., 100., 100., 100., 100., 100.},
/* Z=6 */ {0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 100., 100., 100., 100., 100.},
/* Z=7 */ {0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 100., 100., 100.},
/* Z=8 */ {0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 100.}
  };

  /// \brief Table of chemical element names
  const std::string ParticleTable::elementTable[elementTableSize] = {
    "",
    "H",
    "He",
    "Li",
    "Be",
    "B",
    "C",
    "N",
    "O",
    "F",
    "Ne",
    "Na",
    "Mg",
    "Al",
    "Si",
    "P",
    "S",
    "Cl",
    "Ar",
    "K",
    "Ca",
    "Sc",
    "Ti",
    "V",
    "Cr",
    "Mn",
    "Fe",
    "Co",
    "Ni",
    "Cu",
    "Zn",
    "Ga",
    "Ge",
    "As",
    "Se",
    "Br",
    "Kr",
    "Rb",
    "Sr",
    "Y",
    "Zr",
    "Nb",
    "Mo",
    "Tc",
    "Ru",
    "Rh",
    "Pd",
    "Ag",
    "Cd",
    "In",
    "Sn",
    "Sb",
    "Te",
    "I",
    "Xe",
    "Cs",
    "Ba",
    "La",
    "Ce",
    "Pr",
    "Nd",
    "Pm",
    "Sm",
    "Eu",
    "Gd",
    "Tb",
    "Dy",
    "Ho",
    "Er",
    "Tm",
    "Yb",
    "Lu",
    "Hf",
    "Ta",
    "W",
    "Re",
    "Os",
    "Ir",
    "Pt",
    "Au",
    "Hg",
    "Tl",
    "Pb",
    "Bi",
    "Po",
    "At",
    "Rn",
    "Fr",
    "Ra",
    "Ac",
    "Th",
    "Pa",
    "U",
    "Np",
    "Pu",
    "Am",
    "Cm",
    "Bk",
    "Cf",
    "Es",
    "Fm",
    "Md",
    "No",
    "Lr",
    "Rf",
    "Db",
    "Sg",
    "Bh",
    "Hs",
    "Mt",
    "Ds",
    "Rg",
    "Cn"
  };

  /// \brief Digit names to compose IUPAC element names
  const std::string ParticleTable::elementIUPACDigits = "nubtqphsoe";

  /* Table for cluster decays
   * Definition of "Stable": halflife > 1 ms
   *
   * These table includes decay data for clusters that INCL presently does
   * not produce. It can't hurt.
   *
   * Unphysical nuclides (A<Z) are marked as stable, but should never be
   * produced by INCL. If you find them in the output, something is fishy.
   */
  const ParticleTable::ClusterDecayType ParticleTable::clusterDecayMode[clusterTableZSize][clusterTableASize] =
  {
/*                       A = 0              1               2               3               4                5               6                7               8               9             10             11             12 */
/* Z =  0 */    {StableCluster, StableCluster,   NeutronDecay, NeutronUnbound, NeutronUnbound,  NeutronUnbound, NeutronUnbound,  NeutronUnbound, NeutronUnbound, NeutronUnbound,  NeutronUnbound, NeutronUnbound, NeutronUnbound},
/* Z =  1 */    {StableCluster, StableCluster,  StableCluster,  StableCluster,   NeutronDecay, TwoNeutronDecay,   NeutronDecay, TwoNeutronDecay, NeutronUnbound, NeutronUnbound,  NeutronUnbound, NeutronUnbound, NeutronUnbound},
/* Z =  2 */    {StableCluster, StableCluster,    ProtonDecay,  StableCluster,  StableCluster,    NeutronDecay,  StableCluster,    NeutronDecay,  StableCluster,   NeutronDecay, TwoNeutronDecay, NeutronUnbound, NeutronUnbound},
/* Z =  3 */    {StableCluster, StableCluster,  StableCluster,  ProtonUnbound,    ProtonDecay,     ProtonDecay,  StableCluster,   StableCluster,  StableCluster,  StableCluster,    NeutronDecay,  StableCluster,   NeutronDecay},
/* Z =  4 */    {StableCluster, StableCluster,  StableCluster,  StableCluster,  ProtonUnbound,     ProtonDecay, TwoProtonDecay,   StableCluster,     AlphaDecay,  StableCluster,   StableCluster,  StableCluster,  StableCluster},
/* Z =  5 */    {StableCluster, StableCluster,  StableCluster,  StableCluster,  StableCluster,   ProtonUnbound, TwoProtonDecay,     ProtonDecay,  StableCluster,    ProtonDecay,   StableCluster,  StableCluster,  StableCluster},
/* Z =  6 */    {StableCluster, StableCluster,  StableCluster,  StableCluster,  StableCluster,   StableCluster,  ProtonUnbound,   ProtonUnbound, TwoProtonDecay,  StableCluster,   StableCluster,  StableCluster,  StableCluster},
/* Z =  7 */    {StableCluster, StableCluster,  StableCluster,  StableCluster,  StableCluster,   StableCluster,  StableCluster,   ProtonUnbound,  ProtonUnbound,  ProtonUnbound,     ProtonDecay,    ProtonDecay,  StableCluster},
/* Z =  8 */    {StableCluster, StableCluster,  StableCluster,  StableCluster,  StableCluster,   StableCluster,  StableCluster,   StableCluster,  ProtonUnbound,  ProtonUnbound,   ProtonUnbound,  ProtonUnbound,    ProtonDecay}
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
    ParticleTable::theRealNeutronMass + ParticleTable::theRealChargedPiMass
    + 0.5;
  const G4double ParticleTable::theINCLProtonSeparationEnergy = 6.83;
  const G4double ParticleTable::theINCLNeutronSeparationEnergy = ParticleTable::theINCLProtonSeparationEnergy;
  G4double ParticleTable::protonSeparationEnergy = theINCLProtonSeparationEnergy;
  G4double ParticleTable::neutronSeparationEnergy = theINCLNeutronSeparationEnergy;

#ifdef INCLXX_IN_GEANT4_MODE
  G4IonTable *ParticleTable::theG4IonTable;
#else
  std::vector< std::vector <G4bool> > ParticleTable::massTableMask;
  std::vector< std::vector <G4double> > ParticleTable::massTable;
#endif

  void ParticleTable::initialize(Config const * const theConfig /*=0*/) {
    protonMass = theINCLNucleonMass;
    neutronMass = theINCLNucleonMass;
    piPlusMass = theINCLPionMass;
    piMinusMass = theINCLPionMass;
    piZeroMass = theINCLPionMass;

#ifndef INCLXX_IN_GEANT4_MODE
    std::string dataFilePath;
    if(theConfig)
      dataFilePath = theConfig->getINCLXXDataFilePath();
    readRealMasses(dataFilePath);
#endif

    if(theConfig && theConfig->getUseRealMasses()) {
      getTableMass = ParticleTable::getRealMass;
      getTableParticleMass = ParticleTable::getRealMass;
    } else {
      getTableMass = ParticleTable::getINCLMass;
      getTableParticleMass = ParticleTable::getINCLMass;
    }
#ifdef INCLXX_IN_GEANT4_MODE
    G4ParticleTable *theG4ParticleTable = G4ParticleTable::GetParticleTable();
    theG4IonTable = theG4ParticleTable->GetIonTable();
    theRealProtonMass = theG4ParticleTable->FindParticle("proton")->GetPDGMass() / MeV;
    theRealNeutronMass = theG4ParticleTable->FindParticle("neutron")->GetPDGMass() / MeV;
    theRealChargedPiMass = theG4ParticleTable->FindParticle("pi+")->GetPDGMass() / MeV;
    theRealPiZeroMass = theG4ParticleTable->FindParticle("pi0")->GetPDGMass() / MeV;
#endif

    // Initialise the separation-energy function
    if(!theConfig || theConfig->getSeparationEnergyType()==INCLSeparationEnergy)
      getSeparationEnergy = ParticleTable::getSeparationEnergyINCL;
    else if(theConfig->getSeparationEnergyType()==RealSeparationEnergy)
      getSeparationEnergy = ParticleTable::getSeparationEnergyReal;
    else if(theConfig->getSeparationEnergyType()==RealForLightSeparationEnergy)
      getSeparationEnergy = ParticleTable::getSeparationEnergyRealForLight;
    else {
      FATAL("Unrecognized separation-energy type in ParticleTable initialization: " << theConfig->getSeparationEnergyType() << std::endl);
    }

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

  std::string ParticleTable::getShortName(const ParticleSpecies s) {
    if(s.theType==Composite)
      return getShortName(s.theA,s.theZ);
    else
      return getShortName(s.theType);
  }

  std::string ParticleTable::getName(const ParticleSpecies s) {
    if(s.theType==Composite)
      return getName(s.theA,s.theZ);
    else
      return getName(s.theType);
  }

  std::string ParticleTable::getName(const G4int A, const G4int Z) {
    std::stringstream stream;
    stream << getElementName(Z) << "-" << A;
    return stream.str();
  }

  std::string ParticleTable::getShortName(const G4int A, const G4int Z) {
    std::stringstream stream;
    stream << getElementName(Z);
    if(A>0)
      stream << A;
    return stream.str();
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
    } else if(p == G4INCL::Composite) {
      return std::string("composite");
    }
    return std::string("unknown");
  }

  std::string ParticleTable::getShortName(const ParticleType p) {
    if(p == G4INCL::Proton) {
      return std::string("p");
    } else if(p == G4INCL::Neutron) {
      return std::string("n");
    } else if(p == G4INCL::DeltaPlusPlus) {
      return std::string("d++");
    } else if(p == G4INCL::DeltaPlus) {
      return std::string("d+");
    } else if(p == G4INCL::DeltaZero) {
      return std::string("d0");
    } else if(p == G4INCL::DeltaMinus) {
      return std::string("d-");
    } else if(p == G4INCL::PiPlus) {
      return std::string("pi+");
    } else if(p == G4INCL::PiZero) {
      return std::string("pi0");
    } else if(p == G4INCL::PiMinus) {
      return std::string("pi-");
    } else if(p == G4INCL::Composite) {
      return std::string("comp");
    }
    return std::string("unknown");
  }

  G4double ParticleTable::getINCLMass(const ParticleType pt) {
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

  G4double ParticleTable::getRealMass(const ParticleType t) {
    switch(t) {
      case Proton:
        return theRealProtonMass;
        break;
      case Neutron:
        return theRealNeutronMass;
        break;
      case PiPlus:
      case PiMinus:
        return theRealChargedPiMass;
        break;
      case PiZero:
        return theRealPiZeroMass;
        break;
      default:
        ERROR("Particle::getRealMass : Unknown particle type." << std::endl);
        return 0.0;
        break;
    }
  }

  G4double ParticleTable::getRealMass(const G4int A, const G4int Z) {
// assert(A>=0);
    // For nuclei with Z<0 or Z>A, assume that the exotic charge state is due to pions
    if(Z<0)
      return A*neutronMass - Z*getRealMass(PiMinus);
    else if(Z>A)
      return A*protonMass + (A-Z)*getRealMass(PiPlus);
    else if(Z==0)
      return A*getRealMass(Neutron);
    else if(A==Z)
      return A*getRealMass(Proton);
    else if(A>1) {
#ifndef INCLXX_IN_GEANT4_MODE
      if(ParticleTable::hasMassTable(A,Z))
        return ParticleTable::massTable.at(Z).at(A);
      else {
        DEBUG("Real mass unavailable for isotope A=" << A << ", Z=" << Z
            << ", using Weizsaecker's formula"
            << std::endl);
        return getWeizsaeckerMass(A,Z);
      }
#else
      return theG4IonTable->GetNucleusMass(Z,A) / MeV;
#endif
    } else
      return 0.;
  }

  G4double ParticleTable::getINCLMass(const G4int A, const G4int Z) {
// assert(A>=0);
    // For nuclei with Z<0 or Z>A, assume that the exotic charge state is due to pions
    if(Z<0)
      return A*neutronMass - Z*getINCLMass(PiMinus);
    else if(Z>A)
      return A*protonMass + (A-Z)*getINCLMass(PiPlus);
    else if(A>1)
      return Z*(protonMass - protonSeparationEnergy) + (A-Z)*(neutronMass - neutronSeparationEnergy);
    else if(A==1 && Z==0)
      return getINCLMass(Neutron);
    else if(A==1 && Z==1)
      return getINCLMass(Proton);
    else
      return 0.;
  }

  G4double ParticleTable::getNuclearRadius(const G4int A, const G4int Z) {
// assert(A>0 && Z>=0 && Z<=A);
    if(A >= 19 || (A < 6 && A >= 2)) {
      // For large (Woods-Saxon or Modified Harmonic Oscillator) or small
      // (Gaussian) nuclei, the radius parameter is just the nuclear radius
      return getRadiusParameter(A,Z);
    } else if(A < clusterTableASize && Z < clusterTableZSize && A >= 6) {
      const G4double thisRMS = positionRMS[Z][A];
      if(thisRMS>0.0)
        return thisRMS;
      else {
        ERROR("ParticleTable::getRadiusParameter : Radius for nucleus A = " << A << " Z = " << Z << " is ==0.0" << std::endl);
        return 0.0;
      }
    } else if(A < 19) {
      const G4double theRadiusParameter = getRadiusParameter(A, Z);
      const G4double theDiffusenessParameter = getSurfaceDiffuseness(A, Z);
      // The formula yields the nuclear RMS radius based on the parameters of
      // the nuclear-density function
      return 1.581*theDiffusenessParameter*
        (2.+5.*theRadiusParameter)/(2.+3.*theRadiusParameter);
    } else {
      ERROR("ParticleTable::getNuclearRadius: No radius for nucleus A = " << A << " Z = " << Z << std::endl);
      return 0.0;
    }
  }

  G4double ParticleTable::getRadiusParameter(const G4int A, const G4int Z) {
// assert(A>0 && Z>=0 && Z<=A);
    if(A >= 28) {
      // phenomenological radius fit
      return (2.745e-4 * A + 1.063) * std::pow(A, 1.0/3.0);
    } else if(A < 6 && A >= 2) {
      if(Z<clusterTableZSize) {
        const G4double thisRMS = positionRMS[Z][A];
        if(thisRMS>0.0)
          return thisRMS;
        else {
          ERROR("ParticleTable::getRadiusParameter : Radius for nucleus A = " << A << " Z = " << Z << " is ==0.0" << std::endl);
          return 0.0;
        }
      } else {
        ERROR("ParticleTable::getRadiusParameter : No radius for nucleus A = " << A << " Z = " << Z << std::endl);
        return 0.0;
      }
    } else if(A < 28 && A >= 6) {
      return mediumRadius[A-1];
      //      return 1.581*mediumDiffuseness[A-1]*(2.+5.*mediumRadius[A-1])/(2.+3.*mediumRadius[A-1]);
    } else {
      ERROR("ParticleTable::getRadiusParameter: No radius for nucleus A = " << A << " Z = " << Z << std::endl);
      return 0.0;
    }
  }

  G4double ParticleTable::getMaximumNuclearRadius(const G4int A, const G4int Z) {
    const G4double XFOISA = 8.0;
    if(A >= 19) {
      return getNuclearRadius(A,Z) + XFOISA * getSurfaceDiffuseness(A,Z);
    } else if(A < 19 && A >= 6) {
      return 5.5 + 0.3 * (G4double(A) - 6.0)/12.0;
    } else if(A >= 2) {
      return ParticleTable::getNuclearRadius(A, Z) + 4.5;
    } else {
      ERROR("ParticleTable::getMaximumNuclearRadius : No maximum radius for nucleus A = " << A << " Z = " << Z << std::endl);
      return 0.0;
    }
  }

#ifdef INCLXX_IN_GEANT4_MODE
  G4double ParticleTable::getSurfaceDiffuseness(const G4int A, const G4int /*Z*/ ) {
#else
  G4double ParticleTable::getSurfaceDiffuseness(const G4int A, const G4int Z) {
#endif // INCLXX_IN_GEANT4_MODE

    if(A >= 28) {
      return 1.63e-4 * A + 0.510;
    } else if(A < 28 && A >= 19) {
      return mediumDiffuseness[A-1];
    } else if(A < 19 && A >= 6) {
      return mediumDiffuseness[A-1];
    } else if(A < 6 && A >= 2) {
      ERROR("ParticleTable::getSurfaceDiffuseness: was called for A = " << A << " Z = " << Z << std::endl);
      return 0.0;
    } else {
      ERROR("ParticleTable::getSurfaceDiffuseness: No diffuseness for nucleus A = " << A << " Z = " << Z << std::endl);
      return 0.0;
    }
  }

  std::string ParticleTable::getElementName(const G4int Z) {
    if(Z<1) {
      WARN("getElementName called with Z<1" << std::endl);
      return elementTable[0];
    } else if(Z<elementTableSize)
      return elementTable[Z];
    else
      return getIUPACElementName(Z);
  }

#ifndef INCLXX_IN_GEANT4_MODE
  void ParticleTable::readRealMasses(std::string const &path) {
    // Clear the existing tables, if any
    massTableMask.clear();
    massTable.clear();

    // File name
    std::string fileName(path + "/walletlifetime.dat");
    DEBUG("Reading real nuclear masses from file " << fileName << std::endl);

    // Open the file stream
    std::ifstream massTableIn(fileName.c_str());
    if(!massTableIn.good()) {
      FATAL("Cannot open " << fileName << " data file." << std::endl);
      return;
    }

    // Read in the data
    unsigned int Z, A;
    const G4double amu = 931.494061; // atomic mass unit in MeV/c^2
    const G4double eMass = 0.5109988; // electron mass in MeV/c^2
    G4double excess;
    massTableIn >> A >> Z >> excess;
    do {
      // Dynamically determine the table size
      if(Z>=massTable.size()) {
        massTable.resize(Z+1);
        massTableMask.resize(Z+1);
      }
      if(A>=massTable[Z].size()) {
        massTable[Z].resize(A+1, 0.);
        massTableMask[Z].resize(A+1, false);
      }

      massTable.at(Z).at(A) = A*amu + excess - Z*eMass;
      massTableMask.at(Z).at(A) = true;

      massTableIn >> A >> Z >> excess;
    } while(massTableIn.good());

  }
#endif

  std::string ParticleTable::getIUPACElementName(const G4int Z) {
    std::stringstream elementStream;
    elementStream << Z;
    std::string elementName = elementStream.str();
    std::transform(elementName.begin(), elementName.end(), elementName.begin(), ParticleTable::intToIUPAC);
    elementName[0] = std::toupper(elementName.at(0));
    return elementName;
  }

  G4int ParticleTable::parseIUPACElement(std::string const &s) {
    // Normalise to lower case
    std::string elementName(s);
    std::transform(elementName.begin(), elementName.end(), elementName.begin(), ::tolower);
    // Return 0 if the element name contains anything but IUPAC digits
    if(elementName.find_first_not_of(elementIUPACDigits)!=std::string::npos)
      return 0;
    std::transform(elementName.begin(), elementName.end(), elementName.begin(), ParticleTable::iupacToInt);
    std::stringstream elementStream(elementName);
    G4int Z;
    elementStream >> Z;
    return Z;
  }

}
