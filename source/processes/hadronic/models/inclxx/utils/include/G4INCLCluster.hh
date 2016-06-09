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

#ifndef G4INCLCluster_hh
#define G4INCLCluster_hh 1

#include "G4INCLParticle.hh"

namespace G4INCL {

  /**
   * Cluster is a particle (inherits from the Particle class) that is
   * actually a collection of elementary particles.
   */
  class Cluster : public Particle {
  public:

    /**
     * We can build dummy clusters based on Z and A and set their properties
     * later using standard setters (setEnergy, setPosition...).
     */
    Cluster(const G4int Z, const G4int A) : Particle()
    {
      setType(Composite);
      theZ = Z;
      theA = A;
      setMass(ParticleTable::getMass(theA,theZ));
      setNumberOfCollisions(1);
      makeParticipant();
    };

    /**
     * We can also build cluster starting from a single "leading"
     * particle and then incrementally add more paricles to it using
     * addParticle().
     */
    Cluster(Particle* p) : Particle()
    {
      setType(Composite);
      addParticle(p);
      setMass(ParticleTable::getMass(theA,theZ));
      adjustMomentumFromEnergy();
      setNumberOfCollisions(1);
      makeParticipant();
    };

    /**
     * A cluster can be directly built from a list of particles.
     */
    Cluster(ParticleList *pl) : Particle()
    {
      setType(Composite);
      for(ParticleIter i = pl->begin(); i != pl->end(); ++i) {
	addParticle((*i));
      }
      thePosition /= theA;
      setMass(ParticleTable::getMass(theA,theZ));
      adjustMomentumFromEnergy();
      setNumberOfCollisions(1);
      makeParticipant();
    };

    /**
     * A cluster can be directly built from a list of particles.
     */
    Cluster(const ParticleList &pl) : Particle()
    {
      setType(Composite);
      for(ParticleIter i = pl.begin(); i != pl.end(); ++i) {
	addParticle((*i));
      }
      thePosition /= theA;
      setMass(ParticleTable::getMass(theA,theZ));
      adjustMomentumFromEnergy();
      setNumberOfCollisions(1);
      makeParticipant();
    };

    virtual ~Cluster() {};

    /// \brief Set the charge number of the cluster
    void setZ(const G4int Z) { theZ = Z; }

    /// \brief Set the mass number of the cluster
    void setA(const G4int A) { theA = A; }

    /**
     * Get the list of particles in the cluster.
     */
    ParticleList const *getParticles() const { return &particles; }

    /**
     * Add one particle to the cluster. This updates the cluster mass,
     * energy, size, etc.
     */
    void addParticle(Particle * const p) {
      particles.push_back(p);
      theEnergy += p->getEnergy();
      thePotentialEnergy += p->getPotentialEnergy();
      theMomentum += p->getMomentum();
      thePosition += p->getPosition();
      theA += p->getA();
      theZ += p->getZ();
    };

    /// \brief Returns the list of particles that make up the cluster
    ParticleList getParticleList() const { return particles; }

    std::string prG4int() const {
      std::stringstream ss;
      ss << "Cluster (ID = " << ID << ") type = ";
      ss << ParticleTable::getName(theType);
      ss << std::endl
        << "   energy = " << theEnergy << std::endl
        << "   momentum = "
        << theMomentum.prG4int()
        << std::endl
        << "   position = "
        << thePosition.prG4int()
        << std::endl
        << "Contains the following particles:"
        << std::endl;
      for(ParticleIter i=particles.begin(); i!=particles.end(); ++i)
        ss << (*i)->prG4int();
      return ss.str();
    }

  private:
    ParticleList particles;
  };

}

#endif
