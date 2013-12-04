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

#include "G4INCLStore.hh"
#include <fstream>
#include "G4INCLLogger.hh"
#include "G4INCLCluster.hh"

namespace G4INCL {

  Store::Store(Config const * const config) :
    loadedA(0),
    loadedZ(0),
    loadedStoppingTime(0.),
    theConfig(config)
  {
  }

  Store::~Store() {
    theBook.reset();
    clear();
  }

  void Store::add(Particle *p) {
    inside.push_back(p);

    if(particleAvatarConnections.find(p)==particleAvatarConnections.end()) {
      IAvatarList *avatars = new IAvatarList;
      particleAvatarConnections[p] = avatars;
    }
  }

  void Store::addParticleEntryAvatar(IAvatar *a) {
    // Add the avatar to the avatar map
    avatarList.push_back(a);

    ParticleList pList = a->getParticles();
    for(ParticleIter i=pList.begin(), e=pList.end(); i!=e; ++i) {
      addIncomingParticle((*i));
      // Connect each particle to the avatar
      connectAvatarToParticle(a, *i);
    }
  }

  void Store::addParticleEntryAvatars(IAvatarList const &al) {
    for(IAvatarIter a=al.begin(), e=al.end(); a!=e; ++a)
      addParticleEntryAvatar(*a);
  }

  void Store::add(IAvatar *a) {
    // Add the avatar to the avatar map
    avatarList.push_back(a);

    ParticleList pList = a->getParticles();
    for(ParticleIter i=pList.begin(), e=pList.end(); i!=e; ++i) {
      // If one of the particles participating in this avatar hasn't been
      // registered with the store, it's probably a symptom of a bug
      // somewhere...
// assert(particleAvatarConnections.find(*i) != particleAvatarConnections.end());

      // Connect each particle to the avatar
      connectAvatarToParticle(a, *i);
    }

  }

  void Store::addIncomingParticle(Particle * const p) {
    incoming.push_back(p);
  }

  void Store::connectAvatarToParticle(IAvatar * const a, Particle * const p) {
    std::map<Particle*, IAvatarList*>::const_iterator iter = particleAvatarConnections.find(p);
    // If the particle is already connected to other avatars
    if(iter!=particleAvatarConnections.end()) { // Add to the existing map entry
      iter->second->push_back(a);
    } else { // Create a new map entry
      IAvatarList *avatars = new IAvatarList;
      avatars->push_back(a);
      particleAvatarConnections[p] = avatars;
    }
  }

  void Store::disconnectAvatarFromParticle(IAvatar * const a, Particle * const p) {
    particleAvatarConnections.find(p)->second->remove(a);
  }

  void Store::removeAvatar(IAvatar * const avatar) {
    // Disconnect the avatar from particles
    ParticleList particlesRelatedToAvatar = avatar->getParticles();
    for(ParticleIter particleIter = particlesRelatedToAvatar.begin(), e = particlesRelatedToAvatar.end();
        particleIter != e; ++particleIter) {
      disconnectAvatarFromParticle(avatar, *particleIter);
    }

#ifdef INCL_AVATAR_SEARCH_INCLSort
    // Remove the avatar iterator from the avatarIterList, if it is present.
    std::list<IAvatarIter>::iterator it=binaryIterSearch(avatar);
    if(it != avatarIterList.end())
      avatarIterList.erase(it);
#endif

    // Remove the avatar itself
    avatarList.remove(avatar);
  }

  void Store::removeAndDeleteAvatar(IAvatar * const avatar) {
    removeAvatar(avatar);
    delete avatar;
  }

  void Store::particleHasBeenUpdated(Particle * const particle) {
    // must make a copy of this list, because calls to removeAvatar will modify
    // the list itself
    IAvatarList avatars = *(particleAvatarConnections.find(particle)->second);
    std::for_each(avatars.begin(), avatars.end(), std::bind1st(std::mem_fun(&G4INCL::Store::removeAndDeleteAvatar), this));
  }

#ifdef INCL_AVATAR_SEARCH_INCLSort
  std::list<IAvatarIter>::iterator Store::binaryIterSearch(IAvatar const * const avatar) {
    std::list<IAvatarIter>::iterator it;
    std::iterator_traits<std::list<IAvatarIter>::iterator>::difference_type count, step;
    std::list<IAvatarIter>::iterator first = avatarIterList.begin();
    std::list<IAvatarIter>::iterator last = avatarIterList.end();
    const G4double avatarTime = avatar->getTime();
    count = distance(first,last);
    while (count>0)
    {
      it = first; step=count/2; advance(it,step);
      if ((**it)->getTime()>avatarTime)
      { first=++it; count-=step+1;  }
      else count=step;
    }
    if(first!=last && (**first)->getID()==avatar->getID())
      return first;
    else
      return last;
  }
#endif

  IAvatar* Store::findSmallestTime() {
    if(avatarList.empty()) return NULL;

#ifdef INCL_AVATAR_SEARCH_FullSort

    /* Full sort algorithm.
     *
     * Simple, but guaranteed to work.
     */
    avatarList.sort(Store::avatarComparisonPredicate);
    IAvatar *avatar = avatarList.front();

#elif defined(INCL_AVATAR_SEARCH_INCLSort)

    /* Partial sort algorithm used by INCL4.6.
     *
     * It nevers sorts the whole avatar list, but rather starts from the last
     * best avatar. It requires the avatarList to be updated by appending new
     * avatars at the end.
     */

    IAvatarIter best;
    if(avatarIterList.empty())
      best = avatarList.begin();
    else
      best = avatarIterList.back();
    G4double bestTime = (*best)->getTime();
    IAvatarIter a = best;

    for(++a; a!=avatarList.end(); ++a)
      if((*a)->getTime() < bestTime) {
        best = a;
        bestTime = (*best)->getTime();
        avatarIterList.push_back(best);
      }
    IAvatar *avatar = *best;

#elif defined(INCL_AVATAR_SEARCH_MinElement)

    /* Algorithm provided by the C++ stdlib. */
    IAvatar *avatar = *(std::min_element(avatarList.begin(), avatarList.end(),
          Store::avatarComparisonPredicate));

#else
#error Unrecognized INCL_AVATAR_SEARCH. Allowed values are: FullSort, INCLSort, MinElement.
#endif

    removeAvatar(avatar);
    return avatar;
  }

  void Store::timeStep(G4double step) {
    for(ParticleIter particleIter = inside.begin(), particleEnd=inside.end();
	particleIter != particleEnd; ++particleIter) {
      (*particleIter)->propagate(step);
    }
  }

  void Store::particleHasBeenEjected(Particle * const p) {
    particleHasBeenUpdated(p);
    // The particle will be destroyed when destroying the Store
    inside.remove(p);
    std::map<Particle*, IAvatarList*>::iterator mapItem = particleAvatarConnections.find(p);
    delete mapItem->second;
    particleAvatarConnections.erase(mapItem);
  }

  void Store::particleHasBeenDestroyed(Particle * const p) {
    particleHasBeenUpdated(p);
    // Have to destroy the particle here, the Store will forget about it
    inside.remove(p);
    delete p;
  }

  void Store::particleHasEntered(Particle * const particle) {
    removeFromIncoming(particle);
    add(particle);
  }

  void Store::clearAvatars() {
    for(IAvatarIter iter = avatarList.begin(), e = avatarList.end(); iter != e; ++iter) {
      delete *iter;
    }

    for(std::map<Particle*, IAvatarList*>::iterator iter = particleAvatarConnections.begin(),
	e = particleAvatarConnections.end(); iter != e; ++iter) {
      delete iter->second;
    }

    particleAvatarConnections.clear();
    avatarList.clear();

  }

  void Store::initialiseParticleAvatarConnections() {
    for(ParticleIter ip=inside.begin(), e=inside.end(); ip!=e; ++ip) {
      particleAvatarConnections[*ip] = new IAvatarList;
    }
  }

  void Store::clear() {
    clearAvatars();

    clearInside();
    clearOutgoing();

    if( incoming.size() != 0 ) {
      INCL_WARN("Incoming list is not empty when Store::clear() is called" << std::endl);
    }
    incoming.clear();

#ifdef INCL_AVATAR_SEARCH_INCLSort
    avatarIterList.clear();
#endif

  }

  void Store::clearInside() {
    for(ParticleIter iter=inside.begin(), e=inside.end(); iter!=e; ++iter) {
      delete *iter;
    }
    inside.clear();
  }

  void Store::clearOutgoing() {
    for(ParticleIter iter=outgoing.begin(), e=outgoing.end(); iter!=e; ++iter) {
      if((*iter)->isCluster()) {
        Cluster *c = dynamic_cast<Cluster *>(*iter);
// assert(c);
#ifdef INCLXX_IN_GEANT4_MODE
        if(!c)
          continue;
#endif
        c->deleteParticles();
      }
      delete (*iter);
    }
    outgoing.clear();
  }

  void Store::loadParticles(std::string filename) {
    clear();
    G4int projectileA, projectileZ, A, Z;
    G4double stoppingTime, cutNN;
    G4int ID, type, isParticipant;
    G4double x, y, z;
    G4double px, py, pz, E, v;

    std::ifstream in(filename.c_str());
    in >> projectileA >> projectileZ >> A >> Z >> stoppingTime >> cutNN;
    loadedA = A;
    loadedZ = Z;
    loadedStoppingTime = stoppingTime;

    G4int readA = 0;
    G4int readZ = 0;
    while(1) {
      in >> ID >> type >> isParticipant >> x >> y >> z >> px >> py >> pz >> E >> v;
      if(!in.good()) break;
      ParticleType t;
      if(type == 1) {
	t = Proton;
	readZ++;
	readA++;
      }
      else if(type == -1) {
	t = Neutron;
	readA++;
      }
      else {
        INCL_FATAL("Unrecognized particle type while loading particles; type=" << type << std::endl);
        t = UnknownParticle;
      }

      Particle *p = new Particle(t, E, ThreeVector(px, py, pz),
				 ThreeVector(x, y, z));
      p->setPotentialEnergy(v);
      if(isParticipant == 1) {
        p->makeParticipant();
        theBook.incrementCascading();
      }
      add(p);
    }

    in.close();
  }

  std::string Store::printParticleConfiguration() {
    std::stringstream ss;
    G4int A = 0, Z = 0;
    for(ParticleIter i=inside.begin(), e=inside.end(); i!=e; ++i) {
      if((*i)->getType() == Proton) {
	A++;
	Z++;
      }
      if((*i)->getType() == Neutron) {
	A++;
      }
    }
    // Note: Projectile A and Z are set to 0 (we don't really know
    // anything about them at this point).
    ss << "0 0 " << A << " " << Z << " "
	      << "100.0" << " "
	      << "0.0" << std::endl;

    for(ParticleIter i=inside.begin(), e=inside.end(); i!=e; ++i) {
      G4int ID = (*i)->getID();
      G4int type = 0;
      if((*i)->getType() == Proton) {
	type = 1;
      }
      if((*i)->getType() == Neutron) {
	type = -1;
      }

      G4int isParticipant = 0;
      if((*i)->isParticipant()) {
	isParticipant = 1;
      }

      G4double x = (*i)->getPosition().getX();
      G4double y = (*i)->getPosition().getY();
      G4double z = (*i)->getPosition().getZ();
      G4double E = (*i)->getEnergy();
      G4double px = (*i)->getMomentum().getX();
      G4double py = (*i)->getMomentum().getY();
      G4double pz = (*i)->getMomentum().getZ();
      G4double V = (*i)->getPotentialEnergy();

      ss << ID << " " << type << " " << isParticipant << " "
		<< x << " " << y << " " << z << " "
		<< px << " " << py << " " << pz << " "
	 << E << " " << V << std::endl;
    }

    return ss.str();
  }

  void Store::writeParticles(std::string filename) {
    std::ofstream out(filename.c_str());
    out << printParticleConfiguration();
    out.close();
  }

  std::string Store::printAvatars() {
    std::stringstream ss;
    for(IAvatarIter i = avatarList.begin(), e = avatarList.end(); i != e; ++i) {
      ss << (*i)->toString() << std::endl;
    }
    return ss.str();
  }

  G4bool Store::containsCollisions() const {
    for(IAvatarIter i = avatarList.begin(), e = avatarList.end(); i != e; ++i)
      if((*i)->getType()==CollisionAvatarType) return true;
    return false;
  }

}
