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

#ifndef G4INCLParticleTable_hh
#define G4INCLParticleTable_hh 1

#include <string>
#include <vector>
// #include <cassert>

#include "G4INCLParticleType.hh"
#include "G4INCLParticleSpecies.hh"
#include "G4INCLLogger.hh"
#include "G4INCLConfig.hh"

#ifdef INCLXX_IN_GEANT4_MODE
#include "G4IonTable.hh"
#include "G4ParticleTable.hh"
#endif
#include "G4INCLGlobals.hh"
#include "G4INCLNaturalIsotopicDistributions.hh"

namespace G4INCL {
  class ParticleTable {
  public:
    /// \brief Initialize the particle table
    static void initialize(Config const * const theConfig = 0);

    /// Get the isospin of a particle
    static G4int getIsospin(const ParticleType t);

    /// Get the native INCL name of the particle
    static std::string getName(const ParticleType t);

    /// Get the short INCL name of the particle
    static std::string getShortName(const ParticleType t);

    /// Get the native INCL name of the particle
    static std::string getName(const ParticleSpecies s);

    /// Get the short INCL name of the particle
    static std::string getShortName(const ParticleSpecies s);

    /// Get the native INCL name of the ion
    static std::string getName(const G4int A, const G4int Z);

    /// Get the short INCL name of the ion
    static std::string getShortName(const G4int A, const G4int Z);

    ///\brief Get INCL nuclear mass (in MeV/c^2)
    static G4double getINCLMass(const G4int A, const G4int Z);

    ///\brief Get INCL particle mass (in MeV/c^2)
    static G4double getINCLMass(const ParticleType t);

#ifndef INCLXX_IN_GEANT4_MODE
    ///\brief Do we have this particle mass?
    static G4double hasMassTable(const unsigned int A, const unsigned int Z) {
      return ( Z > 0 && A > 0
          && Z < massTableMask.size() && A < massTableMask.at(Z).size()
          && massTableMask.at(Z).at(A));
    }

    /** \brief Weizsaecker mass formula
     *
     * Return the nuclear mass, as calculated from Weizsaecker's mass formula.
     * Adapted from the Geant4 source.
     *
     * \param A the mass number
     * \param Z the charge number
     * \return the nuclear mass [MeV/c^2]
     */
    static G4double getWeizsaeckerMass(const G4int A, const G4int Z) {
      const G4int Npairing = (A-Z)%2;                  // pairing
      const G4int Zpairing = Z%2;
      const G4double fA = (G4double) A;
      const G4double fZ = (G4double) Z;
      G4double binding =
        - 15.67*fA                          // nuclear volume
        + 17.23*Math::pow23(fA)                // surface energy
        + 93.15*((fA/2.-fZ)*(fA/2.-fZ))/fA       // asymmetry
        + 0.6984523*fZ*fZ*Math::powMinus13(fA);      // coulomb
      if( Npairing == Zpairing ) binding += (Npairing+Zpairing-1) * 12.0 / std::sqrt(fA);  // pairing

      return fZ*getRealMass(Proton)+((G4double)(A-Z))*getRealMass(Neutron)+binding;
    }
#endif

    ///\brief Get particle mass (in MeV/c^2)
    static G4double getRealMass(const G4INCL::ParticleType t);
    ///\brief Get nuclear mass (in MeV/c^2)
    static G4double getRealMass(const G4int A, const G4int Z);

    /**\brief Get Q-value (in MeV/c^2)
     *
     * Uses the getTableMass function to compute the Q-value for the
     * following reaction:
     * \f[ (A_1,Z_1) + (A_2, Z_2) --> (A_1+A_2,Z_1+Z_2) \f]
     */
    static G4double getTableQValue(const G4int A1, const G4int Z1, const G4int A2, const G4int Z2) {
      return getTableMass(A1,Z1) + getTableMass(A2,Z2) - getTableMass(A1+A2,Z1+Z2);
    }

    /**\brief Get Q-value (in MeV/c^2)
     *
     * Uses the getTableMass function to compute the Q-value for the
     * following reaction:
     * \f[ (A_1,Z_1) + (A_2, Z_2) --> (A_3,Z_3) + (A1+A2-A3,Z1+Z2-Z3) \f]
     */
    static G4double getTableQValue(const G4int A1, const G4int Z1, const G4int A2, const G4int Z2, const G4int A3, const G4int Z3) {
      return getTableMass(A1,Z1) + getTableMass(A2,Z2) - getTableMass(A3,Z3) - getTableMass(A1+A2-A3,Z1+Z2-Z3);
    }

    // Typedefs and pointers for transparent handling of mass functions
    typedef G4double (*NuclearMassFn)(const G4int, const G4int);
    typedef G4double (*ParticleMassFn)(const ParticleType);
    static NuclearMassFn getTableMass;
    static ParticleMassFn getTableParticleMass;

    static G4double getTableSpeciesMass(const ParticleSpecies &p) {
      if(p.theType == Composite)
        return (*getTableMass)(p.theA, p.theZ);
      else
        return (*getTableParticleMass)(p.theType);
    }

    // Typedefs and pointers for transparent handling of separation energies
    typedef G4double (*SeparationEnergyFn)(const ParticleType, const G4int, const G4int);
    static SeparationEnergyFn getSeparationEnergy;

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
          return 0;
          break;
      }
    }

    static G4double getNuclearRadius(const G4int A, const G4int Z);
    static G4double getRadiusParameter(const G4int A, const G4int Z);
    static G4double getMaximumNuclearRadius(const G4int A, const G4int Z);
    static G4double getSurfaceDiffuseness(const G4int A, const G4int Z);

    /// \brief Return the RMS of the momentum distribution (light clusters)
    static G4double getMomentumRMS(const G4int A, const G4int Z) {
// assert(Z>=0 && A>=0 && Z<=A);
      if(Z<clusterTableZSize && A<clusterTableASize)
        return momentumRMS[Z][A];
      else
        return Math::sqrtThreeFifths * PhysicalConstants::Pf;
    }

    /// \brief Return INCL's default separation energy
    static G4double getSeparationEnergyINCL(const ParticleType t, const G4int /*A*/, const G4int /*Z*/) {
      if(t==Proton)
        return theINCLProtonSeparationEnergy;
      else if(t==Neutron)
        return theINCLNeutronSeparationEnergy;
      else {
        ERROR("ParticleTable::getSeparationEnergyINCL : Unknown particle type." << std::endl);
        return 0.0;
      }
    }

    /// \brief Return the real separation energy
    static G4double getSeparationEnergyReal(const ParticleType t, const G4int A, const G4int Z) {
      // Real separation energies for all nuclei
      if(t==Proton)
        return (*getTableParticleMass)(Proton) + (*getTableMass)(A-1,Z-1) - (*getTableMass)(A,Z);
      else if(t==Neutron)
        return (*getTableParticleMass)(Neutron) + (*getTableMass)(A-1,Z) - (*getTableMass)(A,Z);
      else {
        ERROR("ParticleTable::getSeparationEnergyReal : Unknown particle type." << std::endl);
        return 0.0;
      }
    }

    /// \brief Return the real separation energy only for light nuclei
    static G4double getSeparationEnergyRealForLight(const ParticleType t, const G4int A, const G4int Z) {
      // Real separation energies for light nuclei, fixed values for heavy nuclei
      if(Z<clusterTableZSize && A<clusterTableASize)
        return getSeparationEnergyReal(t, A, Z);
      else
        return getSeparationEnergyINCL(t, A, Z);
    }

    /// \brief Getter for protonSeparationEnergy
    static G4double getProtonSeparationEnergy() { return protonSeparationEnergy; }

    /// \brief Getter for neutronSeparationEnergy
    static G4double getNeutronSeparationEnergy() { return neutronSeparationEnergy; }

    /// \brief Setter for protonSeparationEnergy
    static void setProtonSeparationEnergy(const G4double s) { protonSeparationEnergy = s; }

    /// \brief Setter for protonSeparationEnergy
    static void setNeutronSeparationEnergy(const G4double s) { neutronSeparationEnergy  = s; }

    /// \brief Get the name of the element from the atomic number
    static std::string getElementName(const G4int Z);
    /// \brief Get the name of an unnamed element from the IUPAC convention
    static std::string getIUPACElementName(const G4int Z);

    /** \brief Parse a IUPAC element name
     *
     * Note: this function is UGLY. Look at it at your own peril.
     *
     * \param pS a normalised string (lowercase)
     * \return the charge number of the nuclide, or zero on fail
     */
    static G4int parseIUPACElement(std::string const &pS);

    const static G4int elementTableSize = 113; // up to Cn

    const static G4double effectiveNucleonMass;
    const static G4double effectiveNucleonMass2;
    const static G4double effectiveDeltaMass;
    const static G4double effectivePionMass;
    const static G4double effectiveDeltaDecayThreshold;

    static const G4int maxClusterMass = 12;
    static const G4int maxClusterCharge = 8;

    const static G4int clusterTableZSize = ParticleTable::maxClusterCharge+1;
    const static G4int clusterTableASize = ParticleTable::maxClusterMass+1;
    const static G4double clusterPosFact[maxClusterMass+1];
    const static G4double clusterPosFact2[maxClusterMass+1];
    const static G4int clusterZMin[maxClusterMass+1]; // Lower limit of Z for cluster of mass A
    const static G4int clusterZMax[maxClusterMass+1]; // Upper limit of Z for cluster of mass A
    const static G4double clusterPhaseSpaceCut[maxClusterMass+1];

#ifdef INCLXX_IN_GEANT4_MODE
    static G4IonTable *theG4IonTable;
#else
    static std::vector< std::vector <G4bool> > massTableMask;
    static std::vector< std::vector <G4double> > massTable;
#endif

    // Enumerator for cluster-decay channels
    enum ClusterDecayType {
      StableCluster,
      NeutronDecay,
      ProtonDecay,
      AlphaDecay,
      TwoProtonDecay,
      TwoNeutronDecay,
      ProtonUnbound,
      NeutronUnbound
    };
    const static ClusterDecayType clusterDecayMode[clusterTableZSize][clusterTableASize];

    /** \brief Coulomb conversion factor, in MeV*fm.
     *
     * \f[ e^2/(4 pi epsilon_0) \f]
     */
    static const G4double eSquared;

    static IsotopicDistribution const &getNaturalIsotopicDistribution(const G4int Z) {
      return getNaturalIsotopicDistributions()->getIsotopicDistribution(Z);
    }

    static G4int drawRandomNaturalIsotope(const G4int Z) {
      return getNaturalIsotopicDistributions()->drawRandomIsotope(Z);
    }

  protected:
    ParticleTable() {};
    ~ParticleTable() {};

  private:
    static const G4double theINCLNucleonMass;
    static const G4double theINCLPionMass;
    static const G4double theINCLNeutronSeparationEnergy;
    static const G4double theINCLProtonSeparationEnergy;
    static G4double protonMass;
    static G4double neutronMass;
    static G4double neutronSeparationEnergy;
    static G4double protonSeparationEnergy;
    static G4double piPlusMass, piMinusMass, piZeroMass;
    static G4double theRealProtonMass;
    static G4double theRealNeutronMass;
    static G4double theRealChargedPiMass;
    static G4double theRealPiZeroMass;

    const static G4int mediumNucleiTableSize = 30;
    const static G4double mediumDiffuseness[mediumNucleiTableSize];
    const static G4double mediumRadius[mediumNucleiTableSize];
    const static G4double positionRMS[clusterTableZSize][clusterTableASize];
    const static G4double momentumRMS[clusterTableZSize][clusterTableASize];

    const static std::string elementTable[elementTableSize];
    
#ifndef INCLXX_IN_GEANT4_MODE
    /// \brief Read nuclear masses from a data file
    static void readRealMasses(std::string const &path);
#endif

    const static std::string elementIUPACDigits;

    /// \brief Transform a IUPAC char to an char representing an integer digit
    static char iupacToInt(char c) {
      return (char)(((G4int)'0')+elementIUPACDigits.find(c));
    }

    /// \brief Transform an integer digit (represented by a char) to a IUPAC char
    static char intToIUPAC(char n) { return elementIUPACDigits.at(n); }

    /// \brief Array of natural isotopic distributions
    static const NaturalIsotopicDistributions *theNaturalIsotopicDistributions;

    /// \brief Get the singleton instance of the natural isotopic distributions
    static const NaturalIsotopicDistributions *getNaturalIsotopicDistributions() {
      if(!theNaturalIsotopicDistributions)
        theNaturalIsotopicDistributions = new NaturalIsotopicDistributions;
      return theNaturalIsotopicDistributions;
    }

  };
}

#endif
