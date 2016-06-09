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

/*
 * Particle.hh
 *
 *  Created on: Jun 5, 2009
 *      Author: Pekka Kaitaniemi
 */

#ifndef PARTICLE_HH_
#define PARTICLE_HH_

#include "G4INCLThreeVector.hh"
#include "G4INCLParticleTable.hh"
#include "G4INCLParticleType.hh"
#include "G4INCLLogger.hh"
#include <list>
#include <sstream>
#include <string>

namespace G4INCL {

  class Particle;

  typedef std::list<G4INCL::Particle*> ParticleList;
  typedef std::list<G4INCL::Particle*>::const_iterator ParticleIter;

  class Particle {
  public:
    Particle();
    Particle(ParticleType t, G4double energy, ThreeVector momentum, ThreeVector position);
    Particle(ParticleType t, ThreeVector momentum, ThreeVector position);
    virtual ~Particle();

    /**
     * Get the particle type.
     * @see G4INCL::ParticleType
     */
    G4INCL::ParticleType getType() const {
      return theType;
    };

    void setType(ParticleType t) {
      theType = t;
      switch(theType)
      {
        case DeltaPlusPlus:
          theA = 1;
          theZ = 2;
          break;
        case Proton:
        case DeltaPlus:
          theA = 1;
          theZ = 1;
          break;
        case Neutron:
        case DeltaZero:
          theA = 1;
          theZ = 0;
          break;
        case DeltaMinus:
          theA = 1;
          theZ = -1;
          break;
        case PiPlus:
          theA = 0;
          theZ = 1;
          break;
        case PiZero:
          theA = 0;
          theZ = 0;
          break;
        case PiMinus:
          theA = 0;
          theZ = -1;
          break;
        case Composite:
         // ERROR("Trying to set particle type to Composite! Construct a Cluster object instead" << std::endl);
          break;
        case UnknownParticle:
          ERROR("Trying to set particle type to Unknown!" << std::endl);
          break;
      }

      if( !isResonance() && t!=Composite )
        theMass = ParticleTable::getMass(theType);
    }

    /**
     * Is this a nucleon?
     */
    G4bool isNucleon() const {
      if(theType == G4INCL::Proton || theType == G4INCL::Neutron)
	return true;
      else
	return false;
    };

    G4bool isParticipant() const {
      return participant;
    }

    void makeParticipant() {
      participant = true;
    }

    void makeSpectator() {
      participant = false;
    }

    /** \brief Is this a pion? */
    G4bool isPion() const { return (theType == PiPlus || theType == PiZero || theType == PiMinus); }

    /** \brief Is it a resonance? */
    inline G4bool isResonance() const { return isDelta(); }

    /** \brief Is it a Delta? */
    inline G4bool isDelta() const {
      return (theType==DeltaPlusPlus || theType==DeltaPlus ||
          theType==DeltaZero || theType==DeltaMinus);
    }

    /** \brief Returns the baryon number. */
    G4int getA() const { return theA; }

    /** \brief Returns the charge number. */
    G4int getZ() const { return theZ; }

    G4double getBeta() const {
      const G4double P = theMomentum.mag();
      return P/theEnergy;
    }

    /**
     * Returns a three vector we can give to the boost() -method.
     *
     * In order to go to the particle rest frame you need to multiply
     * the boost vector by -1.0.
     */
    ThreeVector boostVector() const {
      return theMomentum / theEnergy;
    }

    /**
     * Boost the particle using a boost vector.
     *
     * Example (go to the particle rest frame):
     * particle->boost(particle->boostVector());
     */
    void boost(const ThreeVector &boostVector) {
      const G4double beta2 = boostVector.mag2();
      const G4double gamma = 1.0 / std::sqrt(1.0 - beta2);
      const G4double bp = theMomentum.dot(boostVector);
      const G4double alpha = (gamma*gamma)/(1.0 + gamma);

      theMomentum = theMomentum + boostVector * alpha * bp - boostVector * gamma * theEnergy;
      theEnergy = gamma * (theEnergy - bp);
    }

    /** \brief Get the cached particle mass. */
    inline G4double getMass() const { return theMass; }

    /** \brief Get the the particle invariant mass.
     *
     * Uses the relativistic invariant
     * \f[ m = \sqrt{E^2 - {\vec p}^2}\f]
     **/
    G4double getInvariantMass() const {
      const G4double mass = std::pow(theEnergy, 2) - theMomentum.dot(theMomentum);
      if(mass < 0.0) {
        ERROR("E*E - p*p is negative." << std::endl);
        return 0.0;
      } else {
        return std::sqrt(mass);
      }
    };

    /// \brief Get the particle kinetic energy.
    inline G4double getKineticEnergy() const { return theEnergy - theMass; }

    /// \brief Get the particle potential energy.
    inline G4double getPotentialEnergy() const { return thePotentialEnergy; }

    /// \brief Set the particle potential energy.
    inline void setPotentialEnergy(G4double v) { thePotentialEnergy = v; }

    /**
     * Get the energy of the particle in MeV.
     */
    G4double getEnergy() const
    {
      return theEnergy;
    };

    /**
     * Set the mass of the particle in MeV/c^2.
     */
    void setMass(G4double mass)
    {
      this->theMass = mass;
    }

    /**
     * Set the energy of the particle in MeV.
     */
    void setEnergy(G4double energy)
    {
      this->theEnergy = energy;
    };

    /**
     * Get the momentum vector.
     */
    const G4INCL::ThreeVector &getMomentum() const
    {
      return theMomentum;
    };

    /** Get the angular momentum w.r.t. the origin */
    G4INCL::ThreeVector getAngularMomentum() const
    {
      return thePosition.vector(theMomentum);
    };

    /**
     * Set the momentum vector.
     */
    void setMomentum(const G4INCL::ThreeVector &momentum)
    {
      this->theMomentum = momentum;
    };

    /**
     * Set the position vector.
     */
    const G4INCL::ThreeVector &getPosition() const
    {
      return thePosition;
    };

    void setPosition(const G4INCL::ThreeVector &position)
    {
      this->thePosition = position;
    };

    G4double getHelicity() { return theHelicity; };
    void setHelicity(G4double h) { theHelicity = h; };

    void propagate(G4double step) {
      thePosition += (theMomentum*(step/theEnergy));
    };

    /** \brief Return the number of collisions undergone by the particle. **/
    G4int getNumberOfCollisions() const { return nCollisions; }

    /** \brief Set the number of collisions undergone by the particle. **/
    void setNumberOfCollisions(G4int n) { nCollisions = n; }

    /** \brief Increment the number of collisions undergone by the particle. **/
    void incrementNumberOfCollisions() { nCollisions++; }

    /** \brief Return the number of decays undergone by the particle. **/
    G4int getNumberOfDecays() const { return nDecays; }

    /** \brief Set the number of decays undergone by the particle. **/
    void setNumberOfDecays(G4int n) { nDecays = n; }

    /** \brief Increment the number of decays undergone by the particle. **/
    void incrementNumberOfDecays() { nDecays++; }

    /** \brief Mark the particle as out of its potential well
     *
     * This flag is used to control pions created outside their potential well
     * in delta decay. The pion potential checks it and returns zero if it is
     * true (necessary in order to correctly enforce energy conservation). The
     * Nucleus::applyFinalState() method uses it to determine whether new
     * avatars should be generated for the particle.
     */
    void setOutOfWell() { outOfWell = true; }

    /// \brief Check if the particle is out of its potential well
    G4bool isOutOfWell() const { return outOfWell; }

    void setEmissionTime(G4double t) { emissionTime = t; }
    G4double getEmissionTime() { return emissionTime; };

    /** \brief Transverse component of the position w.r.t. the momentum. */
    ThreeVector getTransversePosition() {
      return thePosition - theMomentum *
        (thePosition.dot(theMomentum)/theMomentum.mag2());
    }

    /** \brief Rescale the momentum to match the total energy. */
    const ThreeVector &adjustMomentumFromEnergy();

    /** \brief Recompute the energy to match the momentum. */
    G4double adjustEnergyFromMomentum();

    /** \brief Check if the particle belongs to a given list **/
    G4bool isInList(ParticleList const &l) const {
      for(ParticleIter i=l.begin(); i!=l.end(); ++i)
        if((*i)->getID()==ID) return true;
      return false;
    }

    G4bool isCluster() const {
      if(theType == Composite) return true;
      else return false;
    }

    std::string prG4int() const {
      std::stringstream ss;
      ss << "Particle (ID = " << ID << ") type = ";
      ss << ParticleTable::getName(theType);
      ss << std::endl
        << "   energy = " << theEnergy << std::endl
        << "   momentum = "
        << theMomentum.prG4int()
        << std::endl
        << "   position = "
        << thePosition.prG4int()
        << std::endl;
      return ss.str();
    };

    std::string dump() const {
      std::stringstream ss;
      ss << "(particle " << ID << " ";
      ss << ParticleTable::getName(theType);
      ss << std::endl
        << thePosition.dump()
        << std::endl
        << theMomentum.dump()
        << std::endl
        << theEnergy << ")" << std::endl;
      return ss.str();
    };

    long getID() const { return ID; };

    /**
     * Return a NULL poG4inter
     */
    ParticleList const *getParticles() const {
      WARN("Particle::getParticles() method was called on a Particle object" << std::endl);
      return 0;
    }

  protected:
    G4int theZ, theA;
    G4bool participant;
    G4INCL::ParticleType theType;
    G4double theEnergy;
    G4INCL::ThreeVector theMomentum;
    G4INCL::ThreeVector thePosition;
    G4int nCollisions;
    G4int nDecays;
    G4double thePotentialEnergy;
    long ID;

  private:
    G4double theHelicity;
    G4double emissionTime;
    G4bool outOfWell;

    G4double theMass;
    static long nextID;



  };
}

#endif /* PARTICLE_HH_ */
