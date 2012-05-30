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
// INCL++ revision: v5.1_rc11
//
#define INCLXX_IN_GEANT4_MODE 1

#include "globals.hh"

#ifndef G4INCLCluster_hh
#define G4INCLCluster_hh 1

#include "G4INCLParticle.hh"
#include "G4INCLNuclearDensity.hh"
#include "G4INCLNuclearDensityFactory.hh"
#include "G4INCLNuclearPotentialConstant.hh"
#include "G4INCLNuclearPotentialIsospin.hh"
#include "G4INCLNuclearPotentialEnergyIsospin.hh"
#include "G4INCLINuclearPotential.hh"
#include "G4INCLConfig.hh"

namespace G4INCL {

  /**
   * Cluster is a particle (inherits from the Particle class) that is
   * actually a collection of elementary particles.
   */
  class Cluster : public Particle {
  public:

    /** \brief Build clusters with a NuclearDensity and a NuclearPotential
     *
     * This constructor should mainly be used when constructing Nucleus or
     * when constructing Clusters to be used as composite projectiles.
     */
    Cluster(const G4int Z, const G4int A, Config const * const conf, const G4bool hardFermiSphere=true) :
      Particle(),
      theDensity(NULL),
      thePotential(NULL),
      theExcitationEnergy(0.),
      theSpin(0.,0.,0.)
    {
      setType(Composite);
      theZ = Z;
      theA = A;
      setINCLMass();

      theDensity = NuclearDensityFactory::createDensity(theA, theZ, hardFermiSphere);

      PotentialType potentialType;
      G4bool pionPotential;
      if(conf) {
        potentialType = conf->getPotentialType();
        pionPotential = conf->getPionPotential();
      } else { // By default we don't use energy dependent
        // potential. This is convenient for some tests.
        potentialType = IsospinPotential;
        pionPotential = true;
      }
      switch(potentialType) {
        case IsospinEnergyPotential:
          thePotential = new NuclearPotential::NuclearPotentialEnergyIsospin(theDensity, pionPotential, hardFermiSphere);
          break;
        case IsospinPotential:
          thePotential = new NuclearPotential::NuclearPotentialIsospin(theDensity, pionPotential, hardFermiSphere);
          break;
        case ConstantPotential:
          thePotential = new NuclearPotential::NuclearPotentialConstant(theDensity, pionPotential, hardFermiSphere);
          break;
        default:
          FATAL("Unrecognized potential type at Nucleus creation." << std::endl);
          std::exit(EXIT_FAILURE);
          break;
      }
    }

    /**
     * We can build dummy clusters based on Z and A and set their properties
     * later using standard setters (setEnergy, setPosition...).
     */
    Cluster(const G4int Z, const G4int A) :
      Particle(),
      theDensity(NULL),
      thePotential(NULL),
      theExcitationEnergy(0.),
      theSpin(0.,0.,0.)
    {
      setType(Composite);
      theZ = Z;
      theA = A;
      setINCLMass();
    }

    /**
     * A cluster can be directly built from a list of particles.
     */
    Cluster(ParticleList *pl) :
      Particle(),
      theDensity(NULL),
      thePotential(NULL),
      theExcitationEnergy(0.),
      theSpin(0.,0.,0.)
    {
      setType(Composite);
      for(ParticleIter i = pl->begin(); i != pl->end(); ++i) {
	addParticle((*i));
      }
      thePosition /= theA;
      setINCLMass();
      adjustMomentumFromEnergy();
    };

    /**
     * A cluster can be directly built from a list of particles.
     */
    Cluster(const ParticleList &pl) :
      Particle(),
      theDensity(NULL),
      thePotential(NULL),
      theExcitationEnergy(0.),
      theSpin(0.,0.,0.)
    {
      setType(Composite);
      for(ParticleIter i = pl.begin(); i != pl.end(); ++i) {
	addParticle((*i));
      }
      thePosition /= theA;
      setINCLMass();
      adjustMomentumFromEnergy();
    };

    virtual ~Cluster() {
      delete thePotential;
    };

    /// \brief Copy constructor
    Cluster(const Cluster &rhs) : Particle(rhs) {
      theDensity = rhs.theDensity;
      thePotential = rhs.thePotential;
      theExcitationEnergy = rhs.theExcitationEnergy;
      deleteParticles();
      for(ParticleIter p=rhs.particles.begin(); p!=rhs.particles.end(); ++p) {
        particles.push_back(new Particle(**p));
      }
    }

    /// \brief Assignment operator
    Cluster &operator=(const Cluster &rhs) {
      Cluster temporaryCluster(rhs);
      Particle::operator=(temporaryCluster);
      swap(temporaryCluster);
      return *this;
    }

    /// \brief Helper method for the assignment operator
    void swap(Cluster &rhs) {
      Particle::swap(rhs);
      std::swap(theExcitationEnergy, rhs.theExcitationEnergy);
      std::swap(theDensity, rhs.theDensity);
      std::swap(thePotential, rhs.thePotential);
      // std::swap is overloaded by std::list and guaranteed to operate in
      // constant time
      std::swap(particles, rhs.particles);
    }

    void deleteParticles() {
      for(ParticleIter p=particles.begin(); p!=particles.end(); ++p) {
        delete (*p);
      }
      clearParticles();
    }

    void clearParticles() { particles.clear(); }

    /// \brief Set the charge number of the cluster
    void setZ(const G4int Z) { theZ = Z; }

    /// \brief Set the mass number of the cluster
    void setA(const G4int A) { theA = A; }

    /// \brief Get the excitation energy of the cluster.
    G4double getExcitationEnergy() const { return theExcitationEnergy; }

    /// \brief Set the excitation energy of the cluster.
    void setExcitationEnergy(const G4double e) { theExcitationEnergy=e; }

    /** \brief Get the real particle mass.
     *
     * Overloads the Particle method.
     */
    inline virtual G4double getTableMass() const { return getRealMass(); }

    /**
     * Get the list of particles in the cluster.
     */
    ParticleList const &getParticles() const { return particles; }

    /// \brief Remove a particle from the cluster components.
    void removeParticle(Particle * const p) { particles.remove(p); }

    /**
     * Add one particle to the cluster. This updates the cluster mass,
     * energy, size, etc.
     */
    void addParticle(Particle * const p) {
// assert(p->isNucleon());
      particles.push_back(p);
      theEnergy += p->getEnergy();
      thePotentialEnergy += p->getPotentialEnergy();
      theMomentum += p->getMomentum();
      thePosition += p->getPosition();
      theA += p->getA();
      theZ += p->getZ();
      nCollisions += p->getNumberOfCollisions();
    };

    /// \brief Returns the list of particles that make up the cluster
    ParticleList getParticleList() const { return particles; }

    std::string print() const {
      std::stringstream ss;
      ss << "Cluster (ID = " << ID << ") type = ";
      ss << ParticleTable::getName(theType);
      ss << std::endl
        << "   A = " << theA << std::endl
        << "   Z = " << theZ << std::endl
        << "   mass = " << getMass() << std::endl
        << "   energy = " << theEnergy << std::endl
        << "   momentum = "
        << theMomentum.print()
        << std::endl
        << "   position = "
        << thePosition.print()
        << std::endl
        << "Contains the following particles:"
        << std::endl;
      for(ParticleIter i=particles.begin(); i!=particles.end(); ++i)
        ss << (*i)->print();
      return ss.str();
    }

    NuclearDensity* getDensity() const { return theDensity; };

    NuclearPotential::INuclearPotential* getPotential() const { return thePotential; };

    /// \brief Initialise the NuclearDensity pointer and sample the particles
    virtual void initializeParticles();

    /// \brief Update the particle potential energy.
    inline void updatePotentialEnergy(Particle *p) {
      p->setPotentialEnergy(thePotential->computePotentialEnergy(p));
    }

    /** \brief Boost to the CM of the component particles
     *
     * The position of all particles in the particles list is shifted so that
     * their centre of mass is in the origin and their total momentum is
     * zero.
     */
    void internalBoostToCM() {

      // First compute the current CM position and total momentum
      ThreeVector theCMPosition, theTotalMomentum;
      G4double theTotalEnergy = 0.0;
      for(ParticleIter p=particles.begin(); p!=particles.end(); ++p) {
        theCMPosition += (*p)->getPosition();
        theTotalMomentum += (*p)->getMomentum();
        theTotalEnergy += (*p)->getEnergy();
      }
      theCMPosition /= theA;
// assert((unsigned int)theA==particles.size());

      // Now determine the CM velocity of the particles
      ThreeVector betaCM = theTotalMomentum / theTotalEnergy;

      // The new particle positions and momenta are scaled by a factor of
      // \f$\sqrt{A/(A-1)}\f$, so that the resulting density distributions in
      // the CM have the same variance as the one we started with.
      const G4double rescaling = std::sqrt(((G4double)theA)/((G4double)(theA-1)));

      // Loop again to boost and reposition
      for(ParticleIter p=particles.begin(); p!=particles.end(); ++p) {
        // We should do the following, but the Fortran version actually does
        // not!
        // (*p)->boost(betaCM);
        // FIXME: here is what the Fortran version does:
        (*p)->setMomentum(((*p)->getMomentum()-theTotalMomentum/theA)*rescaling);

        // Set the CM position of the particles
        (*p)->setPosition(((*p)->getPosition()-theCMPosition)*rescaling);
      }

      // Set the global cluster kinematic variables
      thePosition.setX(0.0);
      thePosition.setY(0.0);
      thePosition.setZ(0.0);
      theMomentum.setX(0.0);
      theMomentum.setY(0.0);
      theMomentum.setZ(0.0);
      theEnergy = getMass();

    }

    /** \brief Put the cluster components off shell
     *
     * The Cluster components are put off shell in such a way that their total
     * energy equals the cluster mass.
     */
    void putParticlesOffShell() {
      // Compute the dynamical potential
      const G4double theDynamicalPotential = computeDynamicalPotential();
      DEBUG("The dynamical potential is " << theDynamicalPotential << " MeV" << std::endl);

      for(ParticleIter p=particles.begin(); p!=particles.end(); ++p) {
        const G4double energy = (*p)->getEnergy() - theDynamicalPotential;
        const ThreeVector &momentum = (*p)->getMomentum();
        // Here particles are put off-shell so that we can satisfy the energy-
        // and momentum-conservation laws
        (*p)->setEnergy(energy);
        (*p)->setMass(std::sqrt(energy*energy - momentum.mag2()));
      }
    }

    /** \brief Set the position of the cluster
     *
     * This overloads the Particle method to take into account that the
     * positions of the cluster members must be updated as well.
     */
    void setPosition(const ThreeVector &position) {
      ThreeVector shift(position-thePosition);
      Particle::setPosition(position);
      for(ParticleIter p=particles.begin(); p!=particles.end(); ++p) {
        (*p)->setPosition((*p)->getPosition()+shift);
      }
    }

    /** \brief Boost the cluster with the indicated velocity
     *
     * The Cluster is boosted as a whole, just like any Particle object;
     * moreover, the internal components (particles list) are also boosted,
     * according to Alain Boudard's off-shell recipe.
     *
     * \param boostVector the velocity to boost to [c]
     */
    void boost(const ThreeVector &boostVector) {
      Particle::boost(boostVector);
      for(ParticleIter p=particles.begin(); p!=particles.end(); ++p) {
        (*p)->boost(boostVector);
        // Apply Lorentz contraction to the particle position
        (*p)->lorentzContract(boostVector,thePosition);
      }
    }

    /** \brief Freeze the internal motion of the particles
     *
     * Each particle is assigned a frozen momentum four-vector determined by
     * the collective cluster velocity. This is used for propagation, but not
     * for dynamics. Normal propagation is restored by calling the
     * Particle::thawPropagation() method, which should be done in
     * InteractionAvatar::postInteraction.
     */
    void freezeInternalMotion() {
      const ThreeVector &normMomentum = theMomentum / getMass();
      for(ParticleIter p=particles.begin(); p!=particles.end(); ++p) {
        const G4double pMass = (*p)->getMass();
        const ThreeVector frozenMomentum = normMomentum * pMass;
        const G4double frozenEnergy = std::sqrt(frozenMomentum.mag2()+pMass*pMass);
        (*p)->setFrozenMomentum(frozenMomentum);
        (*p)->setFrozenEnergy(frozenEnergy);
        (*p)->freezePropagation();
      }
    }

    /** \brief Rotate position and momentum of all the particles
     *
     * This includes the cluster components. Overloads Particle::rotate().
     *
     * \param angle the rotation angle
     * \param axis a unit vector representing the rotation axis
     */
    virtual void rotate(const G4double angle, const ThreeVector &axis) {
      Particle::rotate(angle, axis);
      for(ParticleIter p=particles.begin(); p!=particles.end(); ++p) {
        (*p)->rotate(angle, axis);
      }
    }

    /// \brief Make all the components projectile spectators, too
    virtual void makeProjectileSpectator() {
      Particle::makeProjectileSpectator();
      for(ParticleIter p=particles.begin(); p!=particles.end(); ++p) {
        (*p)->makeProjectileSpectator();
      }
    }

    /// \brief Make all the components target spectators, too
    virtual void makeTargetSpectator() {
      Particle::makeTargetSpectator();
      for(ParticleIter p=particles.begin(); p!=particles.end(); ++p) {
        (*p)->makeTargetSpectator();
      }
    }

    /// \brief Make all the components participants, too
    virtual void makeParticipant() {
      Particle::makeParticipant();
      for(ParticleIter p=particles.begin(); p!=particles.end(); ++p) {
        (*p)->makeParticipant();
      }
    }

    /// \brief Get the spin of the nucleus.
    ThreeVector const &getSpin() const { return theSpin; }

    /// \brief Set the spin of the nucleus.
    void setSpin(const ThreeVector &j) { theSpin = j; }

    /// \brief Get the total angular momentum (orbital + spin)
    G4INCL::ThreeVector getAngularMomentum() const {
      return Particle::getAngularMomentum() + getSpin();
    }

  protected:
    ParticleList particles;
    NuclearDensity *theDensity;
    NuclearPotential::INuclearPotential *thePotential;
    G4double theExcitationEnergy;
    ThreeVector theSpin;

  private:
    /** \brief Compute the dynamical cluster potential
     *
     * Alain Boudard's boost prescription for low-energy beams requires to
     * define a "dynamical potential" that allows us to conserve momentum and
     * energy when boosting the projectile cluster.
     *
     * \param boostVector the velocity to boost to [c]
     */
    G4double computeDynamicalPotential() {
      G4double theDynamicalPotential = 0.0;
      for(ParticleIter p=particles.begin(); p!=particles.end(); ++p) {
        theDynamicalPotential += (*p)->getEnergy();
      }
      theDynamicalPotential -= getTableMass();
      theDynamicalPotential /= theA;

      return theDynamicalPotential;
    }

  };

}

#endif
