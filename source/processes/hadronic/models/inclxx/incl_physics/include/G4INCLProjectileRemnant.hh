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

/** \file G4INCLProjectileRemnant.hh
 * \brief Class for constructing a projectile-like remnant.
 *
 * \date 20 March 2012
 * \author Davide Mancusi
 */

#ifndef G4INCLPROJECTILEREMNANT_HH_
#define G4INCLPROJECTILEREMNANT_HH_

#include "G4INCLCluster.hh"
#include "G4INCLRandom.hh"
#include <vector>
#include <map>
#include <numeric>
#include <functional>

namespace G4INCL {

  /// \brief Helper function for ProjectileRemnant::shuffleStoredComponents
  G4int shuffleComponentsHelper(G4int range);

  class ProjectileRemnant : public Cluster {
    // typedefs for the calculation of the projectile excitation energy
    typedef std::vector<G4double> EnergyLevels;
    typedef std::map<long, G4double> EnergyLevelMap;

    public:
    ProjectileRemnant(ParticleSpecies const species, const G4double kineticEnergy)
      : Cluster(species.theZ, species.theA) {

      // Use the table mass
      setTableMass();

      // Set the kinematics
      const G4double projectileMass = getMass();
      const G4double energy = kineticEnergy + projectileMass;
      const G4double momentumZ = std::sqrt(energy*energy - projectileMass*projectileMass);

      // Initialise the particles
      initializeParticles();
      internalBoostToCM();
      putParticlesOffShell();

      // Store the energy levels of the ProjectileRemnant (used to compute its
      // excitation energy)
      storeEnergyLevels();

      // Boost the whole thing
      const ThreeVector aBoostVector = ThreeVector(0.0, 0.0, momentumZ / energy);
      boost(-aBoostVector);

      // Freeze the internal motion of the particles
      freezeInternalMotion();

      // Set as projectile spectator
      makeProjectileSpectator();
    }

    ~ProjectileRemnant() {
      deleteStoredComponents();
      clearEnergyLevels();
    }

    /// \brief Reset the projectile remnant to the state at the beginning of the cascade
    void reset();

    /** \brief Remove a nucleon from the projectile remnant
     *
     * \param p particle to be removed
     * \param theProjectileCorrection correction to be given to the projectile total energy
     */
    void removeParticle(Particle * const p, const G4double theProjectileCorrection);

    /** \brief Add back dynamical spectators to the projectile remnant
     *
     * Try to add the dynamical spectators back to the projectile remnant.
     * Refuse to do so if this leads to a negative projectile excitation
     * energy.
     *
     * Return a list of rejected dynamical spectators.
     */
    ParticleList addDynamicalSpectators(ParticleList pL);

    /** \brief Add back dynamical spectators to the projectile remnant
     *
     * Try as hard as possible to add back all the dynamical spectators.  Don't
     * add spectators that lead to negative excitation energies. Start by
     * adding all of them, and repeatedly remove the most troublesome one until
     * the excitation energy becomes non-negative.
     *
     * Return a list of rejected dynamical spectators.
     */
    ParticleList addMostDynamicalSpectators(ParticleList pL);

    /// \brief Clear the stored projectile components and delete the particles
    void deleteStoredComponents() {
      for(std::map<long,Particle*>::const_iterator p=storedComponents.begin(); p!=storedComponents.end(); ++p)
        delete p->second;
      clearStoredComponents();
    }

    /// \brief Clear the stored projectile components
    void clearStoredComponents() {
      storedComponents.clear();
    }

    /// \brief Clear the stored energy levels
    void clearEnergyLevels() {
      theInitialEnergyLevels.clear();
      theGroundStateEnergies.clear();
    }

    /** \brief Compute the excitation energy
     *
     * Compute the excitation energy of the projectile-like remnant as the
     * difference between the initial and the present configuration. This
     * follows the algorithm proposed by A. Boudard in INCL4.2-HI, as
     * implemented in Geant4.
     *
     * \return the excitation energy
     */
    G4double computeExcitationEnergy(const long exceptID) const;

    EnergyLevels getPresentEnergyLevels(const long exceptID) const {
      EnergyLevels theEnergyLevels;
      for(ParticleIter p=particles.begin(); p!=particles.end(); ++p) {
        if((*p)->getID()!=exceptID) {
          EnergyLevelMap::const_iterator i = theInitialEnergyLevels.find((*p)->getID());
// assert(i!=theInitialEnergyLevels.end());
          theEnergyLevels.push_back(i->second);
        }
      }
// assert(theEnergyLevels.size()==particles.size()-1);
      return theEnergyLevels;
    }

    /// \brief Store the projectile components
    void storeComponents() {
      for(ParticleIter p=particles.begin(); p!=particles.end(); ++p) {
        // Store the particles (needed for forced CN)
        storedComponents[(*p)->getID()]=new Particle(**p);
      }
    }

    /// \brief Get the number of the stored components
    G4int getNumberStoredComponents() const {
      return storedComponents.size();
    }

    /// \brief Store the energy levels
    void storeEnergyLevels() {
      EnergyLevels energies;

      for(ParticleIter p=particles.begin(); p!=particles.end(); ++p) {
        const G4double theCMEnergy = (*p)->getEnergy();
        // Store the CM energy in the EnergyLevels map
        theInitialEnergyLevels[(*p)->getID()] = theCMEnergy;
        energies.push_back(theCMEnergy);
      }

      std::sort(energies.begin(), energies.end());
// assert(energies.size()==(unsigned int)theA);
      theGroundStateEnergies.resize(energies.size());
      // Compute the partial sums of the CM energies -- they are our reference
      // ground-state energies for any number of nucleons
      std::partial_sum(energies.begin(), energies.end(), theGroundStateEnergies.begin());
    }

    private:

    /// \brief Shuffle the list of stored projectile components
    ParticleList shuffleStoredComponents() {
      ParticleList pL = getStoredComponents();
      std::vector<Particle *> theVector(pL.begin(),pL.end());
      std::random_shuffle(theVector.begin(), theVector.end(), shuffleComponentsHelper);
      return ParticleList(theVector.begin(),theVector.end());
    }

    ParticleList getStoredComponents() const {
      ParticleList pL;
      for(std::map<long,Particle*>::const_iterator p=storedComponents.begin(); p!=storedComponents.end(); ++p)
        pL.push_back(p->second);
      return pL;
    }

    /// \brief Return the stored momentum of a given projectile component
    ThreeVector const &getStoredMomentum(Particle const * const p) const {
      std::map<long,Particle*>::const_iterator i = storedComponents.find(p->getID());
      if(i==storedComponents.end()) {
        ERROR("Couldn't find particle " << p->getID() << " in the list of projectile components" << std::endl);
        return p->getMomentum();
      } else {
        return i->second->getMomentum();
      }
    }

    /** \brief Add back a nucleon to the projectile remnant
     *
     * Try to add a dynamical spectator back to the projectile remnant. Refuse
     * to do so if this leads to a negative projectile excitation energy.
     * Return true on success, false on failure.
     */
    G4bool addDynamicalSpectator(Particle * const p);

    /// \brief Return the stored energy of a given projectile component
    /*    G4double getStoredEnergy(Particle const * const p) {
          std::map<long,Particle*>::const_iterator i = initialProjectileComponents.find(p->getID());
          if(i==initialProjectileComponents.end()) {
          ERROR("Couldn't find particle " << p->getID() << " in the list of projectile components" << std::endl);
          return 0.;
          } else {
          return i->second->getEnergy();
          }
          }*/

    /// \brief Stored projectile components
    std::map<long, Particle*> storedComponents;

    /// \brief Initial energy levels of the projectile
    EnergyLevelMap theInitialEnergyLevels;

    /// \brief Ground-state energies for any number of nucleons
    EnergyLevels theGroundStateEnergies;

  };
}

#endif // G4INCLPROJECTILEREMNANT_HH_

