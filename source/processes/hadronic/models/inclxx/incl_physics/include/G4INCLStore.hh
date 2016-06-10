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

#ifndef G4INCLParticleStore_hh
#define G4INCLParticleStore_hh 1

#include <map>
#include <set>
#include <list>
#include <string>
#include <algorithm>

#include "G4INCLParticle.hh"
#include "G4INCLIAvatar.hh"
#include "G4INCLBook.hh"
#include "G4INCLConfig.hh"

#ifdef INCLXX_IN_GEANT4_MODE
#define INCL_AVATAR_SEARCH_MinElement 1
#endif // INCLXX_IN_GEANT4_MODE

#if !defined(NDEBUG) && !defined(INCLXX_IN_GEANT4_MODE)
// Force instantiation of all the std::multimap<Particle*,IAvatar*> methods for
// debugging purposes
namespace G4INCL {
  class Particle;
  class IAvatar;
}
template class std::multimap<G4INCL::Particle*, G4INCL::IAvatar*>;
#endif

namespace G4INCL {

  /**
   * The purpose of the Store object is to act as a "particle manager"
   * that keeps track ofall the particles in our simulation. It also
   * tracks the avatars and their connections to particles.
   */
  class Store {
  public:
    /**
     * Store constructor
     */
    Store(Config const * const config);

    /**
     * Store destructor
     */
    ~Store();

    /**
     * Add one particle to the store.
     *
     * Particle objects don't know anything about avatars so this
     * method will only do two things:
     * 1. add the particle to the particle map ParticleID -> Particle*
     * 2. add an empty entry for this particle into map AvatarID -> [ParticleID]
     */
    void add(Particle *p);

    /** \brief Add a list of particles to the Store
     *
     * This acts as if add(Particle *) was called on each element of the list.
     */
    void add(ParticleList const &pL);

    /// \brief Add one ParticleEntry avatar
    void addParticleEntryAvatar(IAvatar *a);

    /// \brief Add one ParticleEntry avatar
    void addParticleEntryAvatars(IAvatarList const &al);

    /**
     * Add one avatar to the store
     *
     * Avatars know about the particles they are associated
     * with. Adding an avatar consists of the following steps:
     * 1. Add the new avatar to the avatar list
     * 2. Add any related new particles to the store by calling add(Particle*)
     * (this should not happen, by the time we are adding avatars all particles
     * should have already been added)
     * 3. Connect the particles involved to the avatar in the map:
     * particleAvatarConnections :: ParticleID -> [AvatarID]
     * 4. Add the new avatar to the map:
     * avatarParticleConnections :: AvatarID -> [ParticleID]
     */
    void add(IAvatar *a);

    /**
     * Return the list of avatars
     */
    IAvatarList const &getAvatars() const {
      return avatarList;
    }

    /**
     * Add a particle to the incoming list.
     *
     * \param p particle to add
     */
    void addIncomingParticle(Particle * const p);

    /**
     * Add a particle to the incoming list.
     *
     * \param p particle to add
     */
    void removeFromIncoming(Particle * const p) { incoming.remove(p); }

    /// \brief Clear the incoming list
    inline void clearIncoming() {
      incoming.clear();
    }

    /// \brief Clear the incoming list and delete the particles
    inline void deleteIncoming() {
      for(ParticleIter iter=incoming.begin(), e=incoming.end(); iter!=e; ++iter) {
        delete (*iter);
      }
      clearIncoming();
    }

    /** \brief Notify the Store about a particle update
     *
     * Notify the Store that a particle has been updated. This
     * schedules the removal of obsolete avatars and their disconnection from
     * the particle.
     */
    void particleHasBeenUpdated(Particle * const);

    /// \brief Remove avatars that have been scheduled
    void removeScheduledAvatars();

    /**
     * Find the avatar that has the smallest time.
     */
    IAvatar* findSmallestTime();

    /**
     * Make one time step: propagate particles and subtract the length
     * of the step from the avatar times.
     */
    void timeStep(G4double step);

    /**
     * Mark the particle as ejected. This removes it from the list of
     * inside particles and removes all avatars related to this
     * particle.
     */
    void particleHasBeenEjected(Particle * const);

    /** \brief add the particle to the outgoing particle list.
     *
     * \param p pointer to the particle to be added
     */
    void addToOutgoing(Particle *p) { outgoing.push_back(p); }

    /** \brief Add a list of particles to the outgoing particle list.
     *
     * \param pl list of particles to be added
     */
    void addToOutgoing(ParticleList const &pl) {
      for(ParticleIter p=pl.begin(), e=pl.end(); p!=e; ++p)
        addToOutgoing(*p);
    }

    /**
     * Remove the particle from the system. This also removes all
     * avatars related to this particle.
     */
    void particleHasBeenDestroyed(Particle * const);

    /** \brief Move a particle from incoming to inside
     *
     * \param particle pointer to a particle
     **/
    void particleHasEntered(Particle * const particle);

    /**
     * Return the list of incoming particles (i.e. particles that have yet to
     * enter the cascade).
     */
    ParticleList const & getIncomingParticles() const { return incoming; }

    /**
     * Return the list of outgoing particles (i.e. particles that have left the
     * cascade).
     */
    ParticleList const & getOutgoingParticles() const { return outgoing; }

    /** \brief Returns a list of dynamical spectators
     *
     * Looks in the outgoing list for particles without collisions and decays,
     * removes them from outgoing and returns them in a list.
     *
     * \return the (possibly empty) list of dynamical spectators
     */
    ParticleList extractDynamicalSpectators() {
      ParticleList spectators;
      for(ParticleIter p=outgoing.begin(), e=outgoing.end(); p!=e; ++p) {
        if((*p)->isProjectileSpectator()) {
// assert((*p)->isNucleon());
          spectators.push_back(*p); // add them to the list we will return
        }
      }

      // Now erase them from outgoing
      for(ParticleIter i=spectators.begin(); i!=spectators.end(); ++i) {
        outgoing.remove(*i);
      }

      return spectators;
    }

    /**
     * Return the list of "active" particles (i.e. particles that can
     * participate in collisions).
     */
    ParticleList const & getParticles() const { return inside; }

    /**
     * Return the pointer to the Book object which keeps track of
     * various counters.
     */
    Book &getBook() { return theBook; };

    G4int countCascading() {
      G4int n=0;
      for(ParticleIter i=inside.begin(), e=inside.end(); i!=e; ++i) {
        if(!(*i)->isTargetSpectator())
          ++n;
      }
      return n;
    }

    /**
     * Get the config object
     */
    Config const * getConfig() { return theConfig; };

    /**
     * Clear all avatars and particles from the store.
     *
     * Warning! This actually deletes the objects as well!
     */
    void clear();

    /**
     * Clear all inside particles from the store.
     *
     * Warning! This actually deletes the objects as well!
     */
    void clearInside();

    /**
     * Clear all outgoing particles from the store.
     *
     * Warning! This actually deletes the objects as well!
     */
    void clearOutgoing();

    /**
     * Clear avatars only.
     */
    void clearAvatars();

    /**
     * Load particle configuration from ASCII file (see
     * avatarPredictionTest).
     */
    void loadParticles(std::string const &filename);

    /**
     * Get the value of the nucleus mass number that we read from file
     * with loadParticles.
     */
    G4int getLoadedA() { return loadedA; };

    /**
     * Get the value of the nucleus charge number that we read from file
     * with loadParticles.
     */
    G4int getLoadedZ() { return loadedZ; };

    /**
     * Get the value of the stopping time that we read from file
     * with loadParticles.
     */
    G4double getLoadedStoppingTime() { return loadedStoppingTime; };

    /**
     * Print the nucleon configuration of the nucleus.
     */
    std::string printParticleConfiguration();

    /**
     * Print the nucleon configuration of the nucleus.
     */
    void writeParticles(std::string const &filename);

    /**
     * Print the list of avatars
     */
    std::string printAvatars();

    G4bool containsCollisions() const;

#if defined(INCL_AVATAR_SEARCH_FullSort) || defined(INCL_AVATAR_SEARCH_MinElement)
  /** \brief Comparison predicate for avatars.
   *
   * avatarComparisonPredicate is used by the std::sort or std::min_element
   * functions to compare the avatar objects according to their time.
   *
   * \param lhs pointer to the first avatar
   * \param rhs pointer to the second avatar
   * \return true iff lhs' time is smaller than rhs'.
   */
    static G4bool avatarComparisonPredicate(IAvatar *lhs, IAvatar *rhs) {
      return (lhs->getTime() < rhs->getTime());
    }
#endif

  private:
    /// \brief Dummy copy constructor to shut up Coverity warnings
    Store(const Store &rhs);

    /// \brief Dummy assignment operator to shut up Coverity warnings
    Store &operator=(Store const &rhs);


    /** \brief Connect an avatar to a particle
     *
     * Adds the avatar to the list of avatars where the particle appears. This
     * is typically called when the avatar is created.
     *
     * \param p the particle
     * \param a the avatar
     */
    void connectAvatarToParticle(IAvatar * const a, Particle * const p);

    /** \brief Disconnect an avatar from a particle
     *
     * Removes the avatar from the list of avatars where the particle appears.
     * This is typically called when the avatar has been invalidated or
     * realised.
     *
     * \param p the particle
     * \param a the avatar
     */
    void disconnectAvatarFromParticle(IAvatar * const a, Particle * const p);

    /** \brief Remove an avatar from the list of avatars
     *
     * Removes an avatar from the list of all avatars. The avatar is *not*
     * deleted.
     *
     * \param a the avatar to remove
     */
    void removeAvatar(IAvatar * const a);

  private:
    /**
     * Map particle -> [avatar]
     */
    std::multimap<Particle*, IAvatar*> particleAvatarConnections;
    typedef std::multimap<Particle*, IAvatar*>::value_type PAPair;
    typedef std::multimap<Particle*, IAvatar*>::iterator PAIter;
    typedef std::pair<PAIter, PAIter> PAIterPair;

    /// \brief Set of avatars to be removed
    std::set<IAvatar*> avatarsToBeRemoved;
    typedef std::set<IAvatar*>::const_iterator ASIter;

  private:
    /**
     * List of all avatars
     */
    IAvatarList avatarList;

    /**
     * List of incoming particles
     */
    ParticleList incoming;

    /**
     * List of particles that are inside the nucleus
     */
    ParticleList inside;

    /**
     * List of outgoing particles
     */
    ParticleList outgoing;

    /**
     * List of geometrical spectators
     */
    ParticleList geomSpectators;

    /**
     * The current time in the simulation
     */
    G4double currentTime;

    /**
     * The Book object keeps track of global counters
     */
    Book theBook;

    /**
     * The target nucleus mass number that was loaded from a particle file
     */
    G4int loadedA;

    /**
     * The target nucleus charge number that was loaded from a particle file
     */
    G4int loadedZ;

    /**
     * The stopping time that was loaded from a particle file
     */
    G4double loadedStoppingTime;

    /**
     * Pointer to the Config object
     */
    Config const * theConfig;

  };
}

#endif
