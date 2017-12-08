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

#ifndef G4INCLParticleTable_hh
#define G4INCLParticleTable_hh 1

#include <string>
#include <vector>
// #include <cassert>

#include "G4INCLParticleType.hh"
#include "G4INCLParticleSpecies.hh"
#include "G4INCLLogger.hh"
#include "G4INCLConfig.hh"
#include "G4INCLHFB.hh"

#ifdef INCLXX_IN_GEANT4_MODE
#include "G4IonTable.hh"
#include "G4ParticleTable.hh"
#endif
#include "G4INCLGlobals.hh"
#include "G4INCLNaturalIsotopicDistributions.hh"

namespace G4INCL {

  namespace ParticleTable {

    const G4int maxClusterMass = 12;
    const G4int maxClusterCharge = 8;

    const G4int clusterTableZSize = maxClusterCharge+1;
    const G4int clusterTableASize = maxClusterMass+1;

    const G4double effectiveNucleonMass = 938.2796;
    const G4double effectiveNucleonMass2 = 8.8036860777616e5;
    const G4double effectiveDeltaMass = 1232.0;
    const G4double effectiveDeltaWidth = 130.0;
    const G4double effectivePionMass = 138.0;
    const G4double effectiveLambdaMass = 1115.683;
    const G4double effectiveSigmaMass = 1197.45; // max value
    const G4double effectiveKaonMass = 497.614; // max value
    const G4double effectiveAntiKaonMass = 497.614; // max value
    const G4double effectiveEtaMass = 547.862;
    const G4double effectiveOmegaMass = 782.65;
    const G4double effectiveEtaPrimeMass = 957.78;
    const G4double effectivePhotonMass = 0.0;
    extern G4ThreadLocal G4double minDeltaMass;
    extern G4ThreadLocal G4double minDeltaMass2;
    extern G4ThreadLocal G4double minDeltaMassRndm;

    /// \brief Initialize the particle table
    void initialize(Config const * const theConfig = 0);

    /// \brief Get the isospin of a particle
    G4int getIsospin(const ParticleType t);

    /// \brief Get the native INCL name of the particle
    std::string getName(const ParticleType t);

    /// \brief Get the short INCL name of the particle
    std::string getShortName(const ParticleType t);

    /// \brief Get the native INCL name of the particle
    std::string getName(const ParticleSpecies &s);

    /// \brief Get the short INCL name of the particle
    std::string getShortName(const ParticleSpecies &s);

    /// \brief Get the native INCL name of the ion
    std::string getName(const G4int A, const G4int Z);

    /// \brief Get the short INCL name of the ion
    std::string getShortName(const G4int A, const G4int Z);

    /// \brief Get INCL nuclear mass (in MeV/c^2)
    G4double getINCLMass(const G4int A, const G4int Z);

    /// \brief Get INCL particle mass (in MeV/c^2)
    G4double getINCLMass(const ParticleType t);

#ifndef INCLXX_IN_GEANT4_MODE
    /// \brief Do we have this particle mass?
    G4double hasMassTable(const unsigned int A, const unsigned int Z);

    /** \brief Weizsaecker mass formula
     *
     * Return the nuclear mass, as calculated from Weizsaecker's mass formula.
     * Adapted from the Geant4 source.
     *
     * \param A the mass number
     * \param Z the charge number
     * \return the nuclear mass [MeV/c^2]
     */
    G4double getWeizsaeckerMass(const G4int A, const G4int Z);
#endif

    ///\brief Get particle mass (in MeV/c^2)
    G4double getRealMass(const G4INCL::ParticleType t);
    ///\brief Get nuclear mass (in MeV/c^2)
    G4double getRealMass(const G4int A, const G4int Z);

    /**\brief Get Q-value (in MeV/c^2)
     *
     * Uses the getTableMass function to compute the Q-value for the
     * following reaction:
     * \f[ (A_1,Z_1) + (A_2, Z_2) --> (A_1+A_2,Z_1+Z_2) \f]
     */
    G4double getTableQValue(const G4int A1, const G4int Z1, const G4int A2, const G4int Z2);

    /**\brief Get Q-value (in MeV/c^2)
     *
     * Uses the getTableMass function to compute the Q-value for the
     * following reaction:
     * \f[ (A_1,Z_1) + (A_2, Z_2) --> (A_3,Z_3) + (A1+A2-A3,Z1+Z2-Z3) \f]
     */
    G4double getTableQValue(const G4int A1, const G4int Z1, const G4int A2, const G4int Z2, const G4int A3, const G4int Z3);

    G4double getTableSpeciesMass(const ParticleSpecies &p);

    /// \brief Get mass number from particle type
    G4int getMassNumber(const ParticleType t);

    /// \brief Get charge number from particle type
    G4int getChargeNumber(const ParticleType t);
    
    /// \brief Get strangeness number from particle type
    G4int getStrangenessNumber(const ParticleType t);

    G4double getNuclearRadius(const ParticleType t, const G4int A, const G4int Z);
    G4double getLargestNuclearRadius(const G4int A, const G4int Z);
    G4double getRadiusParameter(const ParticleType t, const G4int A, const G4int Z);
    G4double getMaximumNuclearRadius(const ParticleType t, const G4int A, const G4int Z);
    G4double getSurfaceDiffuseness(const ParticleType t, const G4int A, const G4int Z);

    /// \brief Return the RMS of the momentum distribution (light clusters)
    G4double getMomentumRMS(const G4int A, const G4int Z);

    /// \brief Return INCL's default separation energy
    G4double getSeparationEnergyINCL(const ParticleType t, const G4int /*A*/, const G4int /*Z*/);

    /// \brief Return the real separation energy
    G4double getSeparationEnergyReal(const ParticleType t, const G4int A, const G4int Z);

    /// \brief Return the real separation energy only for light nuclei
    G4double getSeparationEnergyRealForLight(const ParticleType t, const G4int A, const G4int Z);

    /// \brief Getter for protonSeparationEnergy
    G4double getProtonSeparationEnergy();

    /// \brief Getter for neutronSeparationEnergy
    G4double getNeutronSeparationEnergy();

    /// \brief Setter for protonSeparationEnergy
    void setProtonSeparationEnergy(const G4double s);

    /// \brief Setter for protonSeparationEnergy
    void setNeutronSeparationEnergy(const G4double s);

    /// \brief Get the name of the element from the atomic number
    std::string getElementName(const G4int Z);

    /// \brief Get the name of an unnamed element from the IUPAC convention
    std::string getIUPACElementName(const G4int Z);

    /// \brief Get the name of the element from the atomic number
    G4int parseElement(std::string pS);

    /** \brief Parse a IUPAC element name
     *
     * Note: this function is UGLY. Look at it at your own peril.
     *
     * \param pS a normalised string (lowercase)
     * \return the charge number of the nuclide, or zero on fail
     */
    G4int parseIUPACElement(std::string const &pS);

    IsotopicDistribution const &getNaturalIsotopicDistribution(const G4int Z);

    G4int drawRandomNaturalIsotope(const G4int Z);

    // Typedefs and pointers for transparent handling of mass functions
    typedef G4double (*NuclearMassFn)(const G4int, const G4int);
    typedef G4double (*ParticleMassFn)(const ParticleType);
    /// \brief Static pointer to the mass function for nuclei
    extern G4ThreadLocal NuclearMassFn getTableMass;
    /// \brief Static pointer to the mass function for particles
    extern G4ThreadLocal ParticleMassFn getTableParticleMass;

    // Typedefs and pointers for transparent handling of separation energies
    typedef G4double (*SeparationEnergyFn)(const ParticleType, const G4int, const G4int);
    /// \brief Static pointer to the separation-energy function
    extern G4ThreadLocal SeparationEnergyFn getSeparationEnergy;

    // Typedefs and pointers for transparent handling of Fermi momentum
    typedef G4double (*FermiMomentumFn)(const G4int, const G4int);
    extern G4ThreadLocal FermiMomentumFn getFermiMomentum;

    /// \brief Return the constant value of the Fermi momentum
    G4double getFermiMomentumConstant(const G4int /*A*/, const G4int /*Z*/);

    /** \brief Return the constant value of the Fermi momentum - special for light
     *
     * This function should always return PhysicalConstants::Pf for heavy
     * nuclei, and values from the momentumRMS table for light nuclei.
     *
     * \param A mass number
     * \param Z charge number
     */
    G4double getFermiMomentumConstantLight(const G4int A, const G4int Z);

    /** \brief Return the value Fermi momentum from a fit
     *
     * This function returns a fitted Fermi momentum, based on data from Moniz
     * et al., Phys. Rev. Lett. 26 (1971) 445. The fitted functional form is
     * \f[
     * p_F(A)=\alpha-\beta\cdot e^{(-A\cdot\gamma)}
     * \f]
     * with \f$\alpha=259.416\f$ MeV/\f$c\f$, \f$\beta=152.824\f$ MeV/\f$c\f$
     * and \f$\gamma=9.5157\cdot10^{-2}\f$.
     *
     * \param A mass number
     */
    G4double getFermiMomentumMassDependent(const G4int A, const G4int /*Z*/);

    /** \brief Get the value of the r-p correlation coefficient
     *
     * \param t the type of the particle (Proton or Neutron)
     * \return the value of the r-p correlation coefficient
     */
    G4double getRPCorrelationCoefficient(const ParticleType t);

    /// \brief Get the thickness of the neutron skin
    G4double getNeutronSkin();

    /// \brief Get the size of the neutron halo
    G4double getNeutronHalo();

    /// \brief Get the type of pion
    ParticleType getPionType(const G4int isosp);

    /// \brief Get the type of nucleon
    ParticleType getNucleonType(const G4int isosp);

    /// \brief Get the type of delta
    ParticleType getDeltaType(const G4int isosp);

    /// \brief Get the type of sigma
    ParticleType getSigmaType(const G4int isosp);

    /// \brief Get the type of kaon
    ParticleType getKaonType(const G4int isosp);

    /// \brief Get the type of antikaon
    ParticleType getAntiKaonType(const G4int isosp);

    /// \brief Get particle width (in s)
    G4double getWidth(const ParticleType t);
  }
}

#endif

