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
    theBook(new Book),
    loadedA(0),
    loadedZ(0),
    loadedStoppingTime(0.),
    theConfig(config)
  {
  }

  Store::~Store() {
    theBook->reset();
    delete theBook;
    theBook = 0;
    clear();
  }

  void Store::add(Particle *p) {
    const long ID = p->getID();
    particles[ID]=p;
    inside.push_back(p);

    if(particleAvatarConnections.find(ID)==particleAvatarConnections.end()) {
      std::vector<long> *aIDs = new std::vector<long>;
      particleAvatarConnections[ID] = aIDs;
    }
  }

  void Store::addParticleEntryAvatar(IAvatar *a) {
    // Add the avatar to the avatar map
    avatars[a->getID()]=a;
    avatarList.push_back(a);

    ParticleList pList = a->getParticles();
    for(ParticleIter i = pList.begin(); i != pList.end(); ++i) {
      addIncomingParticle((*i));
      // Connect each particle to the avatar
      connectParticleAndAvatar((*i)->getID(), a->getID());
    }
  }

  void Store::addParticleEntryAvatars(IAvatarList const &al) {
    for(IAvatarIter a=al.begin(); a!=al.end(); ++a)
      addParticleEntryAvatar(*a);
  }

  void Store::add(IAvatar *a) {
    // Add the avatar to the avatar map
    avatars[a->getID()]=a;
    avatarList.push_back(a);

    ParticleList pList = a->getParticles();
    for(ParticleIter i = pList.begin(); i != pList.end(); ++i) {
      // If one of the particles participating in this avatar haven't
      // been registered with the store we can do it now. On the other
      // hand, if this happens, it's probably a symptom of a bug
      // somewhere...
      if(particles.find((*i)->getID()) == particles.end()) {
       ERROR("Avatar was added before related particles. This is probably a bug." << std::endl);
        add((*i));
      }
      // Connect each particle to the avatar
      connectParticleAndAvatar((*i)->getID(), a->getID());
    }

  }

  void Store::addIncomingParticle(Particle * const p) {
    incoming.push_back(p);
  }

  void Store::connectParticleAndAvatar(long particleID, long avatarID) {
    std::map<long, std::vector<long>* >::const_iterator iter = particleAvatarConnections.find(particleID);
    // If the particle is already connected to other avatars
    if(iter!=particleAvatarConnections.end()) { // Add to the existing map entry
      std::vector<long> *p = iter->second;
      p->push_back(avatarID);
    } else { // Create a new map entry
      std::vector<long> *p = new std::vector<long>;
      p->push_back(avatarID);
      particleAvatarConnections[particleID]=p;
    }
  }

  void Store::removeAvatarFromParticle(long particleID, long avatarID) {
    std::vector<long>* theseAvatars = particleAvatarConnections.find(particleID)->second;
    std::vector<long>* newAvatars = new std::vector<long>();
    for(std::vector<long>::const_iterator iter = theseAvatars->begin();
	iter != theseAvatars->end(); ++iter) {
      if((*iter) != avatarID) {
	newAvatars->push_back((*iter));
      }
    }
    delete theseAvatars;
    particleAvatarConnections[particleID] = newAvatars;
  }

  void Store::removeAvatarByID(long ID) {
    // Disconnect the avatar from particles
    IAvatar *avatar = avatars.find(ID)->second;
    ParticleList particlesRelatedToAvatar = avatar->getParticles();
    for(ParticleIter particleIDiter
	  = particlesRelatedToAvatar.begin();
	particleIDiter != particlesRelatedToAvatar.end(); ++particleIDiter) {
      removeAvatarFromParticle((*particleIDiter)->getID(), ID);
    }

#ifdef INCL_AVATAR_SEARCH_INCLSort
    // Remove the avatar iterator from the avatarIterList, if it is present.
    std::list<IAvatarIter>::iterator it=binaryIterSearch(avatars.find(ID)->second);
    if(it != avatarIterList.end())
      avatarIterList.erase(it);
#endif

    // Remove the avatar itself
    avatarList.remove(avatar);
    avatars.erase(ID);
  }

  void Store::particleHasBeenUpdated(long particleID) {
    std::vector<long> temp_aIDs;
    std::vector<long> *aIDs = particleAvatarConnections.find(particleID)->second;
    for(std::vector<long>::iterator i = aIDs->begin();
	i != aIDs->end(); ++i) {
      temp_aIDs.push_back((*i));
    }

    for(std::vector<long>::iterator i = temp_aIDs.begin();
	i != temp_aIDs.end(); ++i) {
      IAvatar *tmpAvatar = avatars.find(*i)->second;
      removeAvatarByID((*i));
      delete tmpAvatar;
    }
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

  void Store::removeAvatarsFromParticle(long ID) {
    std::vector<long> *relatedAvatars = particleAvatarConnections.find(ID)->second;
    for(std::vector<long>::const_iterator i = relatedAvatars->begin();
	i != relatedAvatars->end(); ++i) {
      removeAvatarByID((*i));
    }
    delete relatedAvatars;
    particleAvatarConnections.erase(ID);
  }

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

    removeAvatarByID(avatar->getID());
    return avatar;
  }

  void Store::timeStep(G4double step) {
    for(std::map<long, Particle*>::iterator particleIter
	  = particles.begin();
	particleIter != particles.end(); ++particleIter) {
      (*particleIter).second->propagate(step);
    }
  }

  void Store::particleHasBeenEjected(long ID) {
    particleHasBeenUpdated(ID);
    // The particle will be destroyed when destroying the Store
    inside.remove(particles.find(ID)->second);
    particles.erase(ID);
    delete particleAvatarConnections.find(ID)->second;
    particleAvatarConnections.erase(ID);
  }

  void Store::particleHasBeenDestroyed(long ID) {
    particleHasBeenUpdated(ID);
    // Have to destroy the particle here, the Store will forget about it
    Particle * const toDelete = particles.find(ID)->second;
    inside.remove(toDelete);
    delete toDelete;
    particles.erase(ID);
  }

  void Store::particleHasEntered(Particle * const particle) {
    removeFromIncoming(particle);
    add(particle);
  }

  ParticleList Store::getParticipants() {
    WARN("Store::getParticipants is probably slow..." << std::endl);
    ParticleList result;
    for(std::map<long, Particle*>::iterator iter = particles.begin();
	iter != particles.end(); ++iter) {
      if((*iter).second->isParticipant()) {
	result.push_back((*iter).second);
      }
    }
    return result;
  }

  ParticleList Store::getSpectators() {
    WARN("Store::getSpectators is probably slow..." << std::endl);
    ParticleList result;
    for(std::map<long, Particle*>::iterator iter = particles.begin();
	iter != particles.end(); ++iter) {
      if(!(*iter).second->isParticipant()) {
	result.push_back((*iter).second);
      }
    }
    return result;
  }
    
  void Store::clearAvatars() {
    for(std::map<long, IAvatar*>::iterator iter = avatars.begin();
	iter != avatars.end(); ++iter) {
      delete (*iter).second;
    }

    for(std::map<long, std::vector<long>*>::iterator iter = particleAvatarConnections.begin();
	iter != particleAvatarConnections.end(); ++iter) {
      delete (*iter).second;
    }

    particleAvatarConnections.clear();
    avatars.clear();
    avatarList.clear();

  }

  void Store::initialiseParticleAvatarConnections() {
    for(ParticleIter ip=inside.begin(); ip!=inside.end(); ++ip) {
      std::vector<long> *aIDs = new std::vector<long>;
      particleAvatarConnections[(*ip)->getID()] = aIDs;
    }
  }

  void Store::clear() {
    clearAvatars();
    inside.clear();

    for(std::map<long, Particle*>::iterator iter = particles.begin();
	iter != particles.end(); ++iter) {
      delete (*iter).second;
    }
    particles.clear();

    clearOutgoing();

    if( incoming.size() != 0 ) {
      WARN("Incoming list is not empty when Store::clear() is called" << std::endl);
    }
    incoming.clear();

#ifdef INCL_AVATAR_SEARCH_INCLSort
    avatarIterList.clear();
#endif

  }

  void Store::clearOutgoing() {
    for(ParticleIter iter = outgoing.begin();	iter != outgoing.end(); ++iter) {
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
        FATAL("Unrecognized particle type while loading particles; type=" << type << std::endl);
        t = UnknownParticle;
      }

      Particle *p = new Particle(t, E, ThreeVector(px, py, pz),
				 ThreeVector(x, y, z));
      p->setPotentialEnergy(v);
      if(isParticipant == 1) {
        p->makeParticipant();
        theBook->incrementCascading();
      }
      add(p);
    }
    
    in.close();
  }

  std::string Store::printParticleConfiguration() {
    std::stringstream ss;
    G4int A = 0, Z = 0;
    for(ParticleIter i = inside.begin(); i != inside.end(); ++i) {
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

    for(ParticleIter i = inside.begin(); i != inside.end(); ++i) {
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
    IAvatarIter i;
    for(i = avatarList.begin(); i != avatarList.end(); ++i) {
      ss << (*i)->toString() << std::endl;
    }
    return ss.str();
  }

  G4bool Store::containsCollisions() const {
    IAvatarIter i;
    for(i = avatarList.begin(); i != avatarList.end(); ++i)
      if((*i)->getType()==CollisionAvatarType) return true;
    return false;
  }

}
