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

#ifndef G4INCLCluster_hh
#define G4INCLCluster_hh 1

#include "G4INCLParticle.hh"
#include "G4INCLNuclearDensityFactory.hh"
#include "G4INCLParticleSampler.hh"
#include "G4INCLAllocationPool.hh"

namespace G4INCL {

  /**
   * Cluster is a particle (inherits from the Particle class) that is
   * actually a collection of elementary particles.
   */
  class Cluster : public Particle {
  public:

    /** \brief Standard Cluster constructor
     *
     * This constructor should mainly be used when constructing Nucleus or
     * when constructing Clusters to be used as composite projectiles.
     */
    Cluster(const G4int Z, const G4int A, const G4int S, const G4bool createParticleSampler=true) :
      Particle(),
      theExcitationEnergy(0.),
      theSpin(0.,0.,0.),
      theParticleSampler(NULL)
    {
      setType(Composite);
      theZ = Z;
      theA = A;
      theS = S;
      setINCLMass();
      if(createParticleSampler)
        theParticleSampler = new ParticleSampler(A,Z,S);
    }

    /**
     * A cluster can be directly built from a list of particles.
     */
    template<class Iterator>
      Cluster(Iterator begin, Iterator end) :
        Particle(),
        theExcitationEnergy(0.),
        theSpin(0.,0.,0.),
        theParticleSampler(NULL)
    {
      setType(Composite);
      for(Iterator i = begin; i != end; ++i) {
        addParticle(*i);
      }
      thePosition /= theA;
      setINCLMass();
      adjustMomentumFromEnergy();
    }

    virtual ~Cluster() {
      delete theParticleSampler;
    }

    /// \brief Copy constructor
    Cluster(const Cluster &rhs) :
      Particle(rhs),
      theExcitationEnergy(rhs.theExcitationEnergy),
      theSpin(rhs.theSpin)
    {
      for(ParticleIter p=rhs.particles.begin(), e=rhs.particles.end(); p!=e; ++p) {
        particles.push_back(new Particle(**p));
      }
      if(rhs.theParticleSampler)
        theParticleSampler = new ParticleSampler(rhs.theA,rhs.theZ,rhs.theS);
      else
        theParticleSampler = NULL;
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
      std::swap(theSpin, rhs.theSpin);
      // std::swap is overloaded by std::list and guaranteed to operate in
      // constant time
      std::swap(particles, rhs.particles);
      std::swap(theParticleSampler, rhs.theParticleSampler);
    }

    ParticleSpecies getSpecies() const {
      return ParticleSpecies(theA, theZ, theS);
    }

    void deleteParticles() {
      for(ParticleIter p=particles.begin(), e=particles.end(); p!=e; ++p) {
        delete (*p);
      }
      clearParticles();
    }

    void clearParticles() { particles.clear(); }

    /// \brief Set the charge number of the cluster
    void setZ(const G4int Z) { theZ = Z; }

    /// \brief Set the mass number of the cluster
    void setA(const G4int A) { theA = A; }

    /// \brief Set the strangess number of the cluster
    void setS(const G4int S) { theS = S; }

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
      particles.push_back(p);
      theEnergy += p->getEnergy();
      thePotentialEnergy += p->getPotentialEnergy();
      theMomentum += p->getMomentum();
      thePosition += p->getPosition();
      theA += p->getA();
      theZ += p->getZ();
      theS += p->getS();
      nCollisions += p->getNumberOfCollisions();
    }

    /// \brief Set total cluster mass, energy, size, etc. from the particles
    void updateClusterParameters() {
      theEnergy = 0.;
      thePotentialEnergy = 0.;
      theMomentum = ThreeVector();
      thePosition = ThreeVector();
      theA = 0;
      theZ = 0;
      theS = 0;
      nCollisions = 0;
      for(ParticleIter p=particles.begin(), e=particles.end(); p!=e; ++p) {
        theEnergy += (*p)->getEnergy();
        thePotentialEnergy += (*p)->getPotentialEnergy();
        theMomentum += (*p)->getMomentum();
        thePosition += (*p)->getPosition();
        theA += (*p)->getA();
        theZ += (*p)->getZ();
        theS += (*p)->getS();
        nCollisions += (*p)->getNumberOfCollisions();
      }
    }

    /// \brief Add a list of particles to the cluster
    void addParticles(ParticleList const &pL) {
      particles = pL;
      updateClusterParameters();
    }

    /// \brief Returns the list of particles that make up the cluster
    ParticleList getParticleList() const { return particles; }

    std::string print() const {
      std::stringstream ss;
      ss << "Cluster (ID = " << ID << ") type = ";
      ss << ParticleTable::getName(theType);
      ss << '\n'
        << "   A = " << theA << '\n'
        << "   Z = " << theZ << '\n'
        << "   S = " << theS << '\n'
        << "   mass = " << getMass() << '\n'
        << "   energy = " << theEnergy << '\n'
        << "   momentum = "
        << theMomentum.print()
        << '\n'
        << "   position = "
        << thePosition.print()
        << '\n'
        << "Contains the following particles:"
        << '\n';
      for(ParticleIter i=particles.begin(), e=particles.end(); i!=e; ++i)
        ss << (*i)->print();
      ss << '\n';
      return ss.str();
    }

    /// \brief Initialise the NuclearDensity pointer and sample the particles
    virtual void initializeParticles();

    /** \brief Boost to the CM of the component particles
     *
     * The position of all particles in the particles list is shifted so that
     * their centre of mass is in the origin and their total momentum is
     * zero.
     */
    void internalBoostToCM() {

      // First compute the current CM position and total momentum
      ThreeVector theCMPosition, theTotalMomentum;
//    G4double theTotalEnergy = 0.0;
      for(ParticleIter p=particles.begin(), e=particles.end(); p!=e; ++p) {
        theCMPosition += (*p)->getPosition();
        theTotalMomentum += (*p)->getMomentum();
//      theTotalEnergy += (*p)->getEnergy();
      }
      theCMPosition /= theA;
// assert((unsigned int)theA==particles.size());

      // Now determine the CM velocity of the particles
      // commented out because currently unused, see below
      // ThreeVector betaCM = theTotalMomentum / theTotalEnergy;

      // The new particle positions and momenta are scaled by a factor of
      // \f$\sqrt{A/(A-1)}\f$, so that the resulting density distributions in
      // the CM have the same variance as the one we started with.
      const G4double rescaling = std::sqrt(((G4double)theA)/((G4double)(theA-1)));

      // Loop again to boost and reposition
      for(ParticleIter p=particles.begin(), e=particles.end(); p!=e; ++p) {
        // \bug{We should do the following, but the Fortran version actually
        // does not!
        // (*p)->boost(betaCM);
        // Here is what the Fortran version does:}
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

      INCL_DEBUG("Cluster boosted to internal CM:" << '\n' << print());

    }

    /** \brief Put the cluster components off shell
     *
     * The Cluster components are put off shell in such a way that their total
     * energy equals the cluster mass.
     */
    void putParticlesOffShell() {
      // Compute the dynamical potential
      const G4double theDynamicalPotential = computeDynamicalPotential();
      INCL_DEBUG("The dynamical potential is " << theDynamicalPotential << " MeV" << '\n');

      for(ParticleIter p=particles.begin(), e=particles.end(); p!=e; ++p) {
        const G4double energy = (*p)->getEnergy() - theDynamicalPotential;
        const ThreeVector &momentum = (*p)->getMomentum();
        // Here particles are put off-shell so that we can satisfy the energy-
        // and momentum-conservation laws
        (*p)->setEnergy(energy);
        (*p)->setMass(std::sqrt(energy*energy - momentum.mag2()));
      }
      INCL_DEBUG("Cluster components are now off shell:" << '\n'
            << print());
    }

    /** \brief Set the position of the cluster
     *
     * This overloads the Particle method to take into account that the
     * positions of the cluster members must be updated as well.
     */
    void setPosition(const ThreeVector &position) {
      ThreeVector shift(position-thePosition);
      Particle::setPosition(position);
      for(ParticleIter p=particles.begin(), e=particles.end(); p!=e; ++p) {
        (*p)->setPosition((*p)->getPosition()+shift);
      }
    }

    /** \brief Boost the cluster with the indicated velocity
     *
     * The Cluster is boosted as a whole, just like any Particle object;
     * moreover, the internal components (particles list) are also boosted,
     * according to Alain Boudard's off-shell recipe.
     *
     * \param aBoostVector the velocity to boost to [c]
     */
    void boost(const ThreeVector &aBoostVector) {
      Particle::boost(aBoostVector);
      for(ParticleIter p=particles.begin(), e=particles.end(); p!=e; ++p) {
        (*p)->boost(aBoostVector);
        // Apply Lorentz contraction to the particle position
        (*p)->lorentzContract(aBoostVector,thePosition);
        (*p)->rpCorrelate();
      }

      INCL_DEBUG("Cluster was boosted with (bx,by,bz)=("
          << aBoostVector.getX() << ", " << aBoostVector.getY() << ", " << aBoostVector.getZ() << "):"
          << '\n' << print());

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
      for(ParticleIter p=particles.begin(), e=particles.end(); p!=e; ++p) {
        const G4double pMass = (*p)->getMass();
        const ThreeVector frozenMomentum = normMomentum * pMass;
        const G4double frozenEnergy = std::sqrt(frozenMomentum.mag2()+pMass*pMass);
        (*p)->setFrozenMomentum(frozenMomentum);
        (*p)->setFrozenEnergy(frozenEnergy);
        (*p)->freezePropagation();
      }
    }

    /** \brief Rotate position of all the particles
     *
     * This includes the cluster components. Overloads Particle::rotateMomentum().
     *
     * \param angle the rotation angle
     * \param axis a unit vector representing the rotation axis
     */
    virtual void rotatePosition(const G4double angle, const ThreeVector &axis);

    /** \brief Rotate momentum of all the particles
     *
     * This includes the cluster components. Overloads Particle::rotateMomentum().
     *
     * \param angle the rotation angle
     * \param axis a unit vector representing the rotation axis
     */
    virtual void rotateMomentum(const G4double angle, const ThreeVector &axis);

    /// \brief Make all the components projectile spectators, too
    virtual void makeProjectileSpectator() {
      Particle::makeProjectileSpectator();
      for(ParticleIter p=particles.begin(), e=particles.end(); p!=e; ++p) {
        (*p)->makeProjectileSpectator();
      }
    }

    /// \brief Make all the components target spectators, too
    virtual void makeTargetSpectator() {
      Particle::makeTargetSpectator();
      for(ParticleIter p=particles.begin(), e=particles.end(); p!=e; ++p) {
        (*p)->makeTargetSpectator();
      }
    }

    /// \brief Make all the components participants, too
    virtual void makeParticipant() {
      Particle::makeParticipant();
      for(ParticleIter p=particles.begin(), e=particles.end(); p!=e; ++p) {
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

  private:
    /** \brief Compute the dynamical cluster potential
     *
     * Alain Boudard's boost prescription for low-energy beams requires to
     * define a "dynamical potential" that allows us to conserve momentum and
     * energy when boosting the projectile cluster.
     */
    G4double computeDynamicalPotential() {
      G4double theDynamicalPotential = 0.0;
      for(ParticleIter p=particles.begin(), e=particles.end(); p!=e; ++p) {
        theDynamicalPotential += (*p)->getEnergy();
      }
      theDynamicalPotential -= getTableMass();
      theDynamicalPotential /= theA;

      return theDynamicalPotential;
    }

  protected:
    ParticleList particles;
    G4double theExcitationEnergy;
    ThreeVector theSpin;
    ParticleSampler *theParticleSampler;

    INCL_DECLARE_ALLOCATION_POOL(Cluster)
  };

}

#endif
