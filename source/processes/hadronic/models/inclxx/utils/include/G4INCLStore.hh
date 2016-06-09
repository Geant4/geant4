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

#ifndef G4INCLParticleStore_hh
#define G4INCLParticleStore_hh 1

#include <map>
#include <list>
#include <vector>
#include <string>
#include <algorithm>

#include "G4INCLParticle.hh"
#include "G4INCLIAvatar.hh"
#include "G4INCLBook.hh"
#include "G4INCLConfig.hh"

#ifdef INCLXX_IN_GEANT4_MODE
#define INCL_AVATAR_SEARCH_MinElement 1
#endif // INCLXX_IN_GEANT4_MODE

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
     * 2. add an empty entry for this particle G4into map AvatarID -> [ParticleID]
     */
    void add(Particle *p);

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
    std::list<IAvatar*> getAvatars() const {
      return avatarList;
    }

    /**
     * Add a particle to the incoming list.
     *
     * \param p particle to add
     */
    void addIncomingParticle(Particle * const p);

    /**
     * Notify the Store that a particle has been updated. This
     * triggers the removal of obsolete avatars and their
     * disconnection from the particle.
     */
    void particleHasBeenUpdated(long);

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
    void particleHasBeenEjected(long);

    /** \brief Add the particle to the outgoing particle list.
     *
     * \param p poG4inter to the particle to be added
     */
    void addToOutgoing(Particle *p) { outgoing.push_back(p); }

    /**
     * Remove the particle from the system. This also removes all
     * avatars related to this particle.
     */
    void particleHasBeenDestroyed(long);

    /** \brief Move a particle from incoming to inside
     *
     * \param particle poG4inter to a particle
     **/
    void particleHasEntered(Particle * const particle);

    /** \brief Return the number of incoming particles.
     *
     * WARNING: this is not the running size of the incoming list! The incoming
     * list is supposed to be filled once and for all, before the cascade. The
     * value returned by this method is the size of the list at that time, i.e.
     * the total number of incoming particles.
     */
    G4int getNumberOfIncomingParticles() const { return nIncomingParticles; }

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

    /**
     * Return the list of "active" particles (i.e. particles that can
     * participate in collisions).
     */
    ParticleList const & getParticles() const { return inside; }

    /**
     * Return the poG4inter to the Book object which keeps track of
     * various counters.
     */
    Book* getBook() {return theBook; };

    G4int countParticipants() {
      G4int n=0;
      for(ParticleIter i=inside.begin(); i!=inside.end(); ++i) {
        if((*i)->isParticipant())
          ++n;
      }
      return n;
    }

    /**
     * Set the config object
     */
    //    void setConfig(Config const * const c) { theConfig = c; };

    /**
     * Get the config object
     */
    Config const * getConfig() { return theConfig; };

    /**
     * Get list of participants (active nucleons).
     *
     * Warning: This (slow) method may be deprecated in the near future...
     */
    ParticleList getParticipants();

    /**
     * Get list of spectators (active nucleons).
     *
     * Warning: This (slow) method may be deprecated in the near future...
     */
    ParticleList getSpectators();

    /**
     * Clear all avatars and particles from the store.
     *
     * Warning! This actually deletes the objects as well!
     */
    void clear();

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

    /** \brief Initialise the particleAvatarConnections map
     *
     * Generate an empty avatar-ID vector for each particle in the inside list
     * and fill in the relevant particle-avatar map entry.
     */
    void initialiseParticleAvatarConnections();

    /**
     * Load particle configuration from ASCII file (see
     * avatarPredictionTest).
     */
    void loadParticles(std::string filename);

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
     * PrG4int the nucleon configuration of the nucleus.
     */
    std::string prG4intParticleConfiguration();

    /**
     * PrG4int the nucleon configuration of the nucleus.
     */
    void writeParticles(std::string filename);

    /**
     * PrG4int the list of avatars
     */
    std::string prG4intAvatars();

    G4bool containsCollisions() const;

#if defined(INCL_AVATAR_SEARCH_FullSort) || defined(INCL_AVATAR_SEARCH_MinElement)
  /** \brief Comparison predicate for avatars.
   *
   * avatarComparisonPredicate is used by the std::sort or std::min_element
   * functions to compare the avatar objects according to their time.
   *
   * \param lhs poG4inter to the first avatar
   * \param rhs poG4inter to the second avatar
   * \return true iff lhs' time is smaller than rhs'.
   */
    static G4bool avatarComparisonPredicate(IAvatar *lhs, IAvatar *rhs) {
      return (lhs->getTime() < rhs->getTime());
    }
#elif defined(INCL_AVATAR_SEARCH_INCLSort)
    /** \brief Perform a binary search on the avatarIterList.
     *
     * By construction, the avatarIterList is always sorted in descending time
     * order. Thus, we can use binary search if we need to look for a specific
     * avatar in the list.
     *
     * Adapted from  STL's binary_search algorithm, as seen on
     * http://www.cplusplus.com/reference/algorithm/binary_search/.
     *
     * \param avatar a poG4inter to the searched avatar.
     * \return an iterator to the IAvatarIter, if the avatar is found; otherwise,
     *         IAvatarList.end().
     */
    std::list<IAvatarIter>::iterator binaryIterSearch(IAvatar const * const avatar);
#endif

  private:
    /**
     * Remove all avatars connedted to a particle
     */
    void removeAvatarsFromParticle(long ID);

    /**
     * Check if a particle is in the particleAvatarConnections map
     */
    G4bool particleInConnectionMap(long);

    /**
     * Check if a particle is in the avatarParticleConnections map
     */
    G4bool avatarInConnectionMap(long);


    /**
     * Connects a particle and an avatar
     */
    void connectParticleAndAvatar(long particleID, long avatarID);

    /**
     * Disconnects a particle and an avatar
     */
    void disconnectAvatarFromParticle(long ID);

    /**
     * Removes an avatar 
     */
    void removeAvatarFromParticle(long particleID, long avatarID);
    void removeAvatarByID(long ID);

  private:
    /**
     * Map of particle ID -> Particle*
     */
    std::map<long, Particle*> particles;

    /**
     * Map of avatar ID -> IAvatar*
     */
    std::map<long, IAvatar*> avatars;

    /**
     * Map particle ID -> [avatar IDs]
     */
    std::map<long, std::vector<long>* > particleAvatarConnections;

    /**
     * List of all avatars
     */
    std::list<IAvatar*> avatarList;

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
     * The current time in the simulation
     */
    G4double currentTime;

    /**
     * The number of incoming particles
     */
    G4int nIncomingParticles;

    /**
     * The Book object keeps track of global counters
     */
    Book *theBook;

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
     * PoG4inter to the Config object
     */
    Config const * const theConfig;

#ifdef INCL_AVATAR_SEARCH_INCLSort
    /** \brief Internal stack for the INCLSort search algorithm.
     *
     * List of std::list<IAvatar*>::const_iterator to keep track of the best
     * avatars so far.
     */
    std::list<IAvatarIter> avatarIterList;
#endif
  };
}
#endif
