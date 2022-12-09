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
  }

  void Store::add(ParticleList const &pL) {
    inside.insert(inside.end(), pL.begin(), pL.end());
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
      // Connect each particle to the avatar
      connectAvatarToParticle(a, *i);
    }

  }

  void Store::addIncomingParticle(Particle * const p) {
    incoming.push_back(p);
  }

  void Store::connectAvatarToParticle(IAvatar * const a, Particle * const p) {
    particleAvatarConnections.insert(PAPair(p,a));
  }

  void Store::disconnectAvatarFromParticle(IAvatar * const a, Particle * const p) {
    PAIterPair iterPair = particleAvatarConnections.equal_range(p);
    for(PAIter i=iterPair.first, last=iterPair.second; i!=last; ++i) {
      if(i->second==a) {
        particleAvatarConnections.erase(i);
        return;
      }
    }
    INCL_WARN("Loop in Store::disconnectAvatarFromParticle fell through." << std::endl
              << "This indicates an inconsistent state of the particleAvatarConnections map." << std::endl);
  }

  void Store::removeAvatar(IAvatar * const avatar) {
    // Disconnect the avatar from particles
    ParticleList particlesRelatedToAvatar = avatar->getParticles();
    for(ParticleIter particleIter = particlesRelatedToAvatar.begin(), e = particlesRelatedToAvatar.end();
        particleIter != e; ++particleIter) {
      disconnectAvatarFromParticle(avatar, *particleIter);
    }

    // Remove the avatar itself
    avatarList.remove(avatar);
  }

  void Store::particleHasBeenUpdated(Particle * const particle) {
    PAIterPair iterPair = particleAvatarConnections.equal_range(particle);
    for(PAIter i=iterPair.first, last=iterPair.second; i!=last; ++i) {
      avatarsToBeRemoved.insert(i->second);
    }
  }

  void Store::removeScheduledAvatars() {
    for(ASIter a=avatarsToBeRemoved.begin(), e=avatarsToBeRemoved.end(); a!=e; ++a) {
      removeAvatar(*a);
      delete *a;
    }
    avatarsToBeRemoved.clear();
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

#elif defined(INCL_AVATAR_SEARCH_MinElement)

    /* Algorithm provided by the C++ stdlib. */
    IAvatar *avatar = *(std::min_element(avatarList.begin(), avatarList.end(),
          Store::avatarComparisonPredicate));

#else
#error Unrecognized INCL_AVATAR_SEARCH. Allowed values are: FullSort, MinElement.
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

    particleAvatarConnections.clear();
    avatarList.clear();
    avatarsToBeRemoved.clear();
  }

  void Store::clear() {
    clearAvatars();

    clearInside();
    clearOutgoing();

    if( incoming.size() != 0 ) {
      INCL_WARN("Incoming list is not empty when Store::clear() is called" << '\n');
    }
    incoming.clear();

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

  void Store::loadParticles(std::string const &filename) {
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

    while(1) { /* Loop checking, 10.07.2015, D.Mancusi */
      in >> ID >> type >> isParticipant >> x >> y >> z >> px >> py >> pz >> E >> v;
      if(!in.good()) break;
      ParticleType t;
      if(type == 1) {
	t = Proton;
      }
      else if(type == -1) {
	t = Neutron;
      }
      else {
        INCL_FATAL("Unrecognized particle type while loading particles; type=" << type << '\n');
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
	      << "0.0" << '\n';

    for(ParticleIter i=inside.begin(), e=inside.end(); i!=e; ++i) {
      G4int ID = (G4int)(*i)->getID();
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
	 << E << " " << V << '\n';
    }

    return ss.str();
  }

  void Store::writeParticles(std::string const &filename) {
    std::ofstream out(filename.c_str());
    out << printParticleConfiguration();
    out.close();
  }

  std::string Store::printAvatars() {
    std::stringstream ss;
    for(IAvatarIter i = avatarList.begin(), e = avatarList.end(); i != e; ++i) {
      ss << (*i)->toString() << '\n';
    }
    return ss.str();
  }

  G4bool Store::containsCollisions() const {
    for(IAvatarIter i = avatarList.begin(), e = avatarList.end(); i != e; ++i)
      if((*i)->getType()==CollisionAvatarType) return true;
    return false;
  }

}
