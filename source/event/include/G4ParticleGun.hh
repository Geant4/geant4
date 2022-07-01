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
// G4ParticleGum
//
// Class description:
//
// This is a concrete class of G4VPrimaryGenerator. It shoots a particle of
// given type into a given direction with either a given kinetic energy or
// momentum.
// The position and time of the primary particle must be set by the
// corresponding set methods of the G4VPrimaryGenerator base class, otherwise
// zero will be set.
//
// The FAQ to this class is for randomizing position/direction/kinetic energy
// of the primary particle. But, G4ParticleGun does NOT have any way of
// randomization. Instead, the user's concrete implementation of
// G4VUserPrimaryGeneratorAction which transfer the G4Event object
// to this particle gun can randomize these quantities and set to this
// particle gun before invoking GeneratePrimaryVertex() method.
// Note that, even if the particle gun shoots more than one particles at one
// invokation of the GeneratePrimaryVertex() method, all particles have the
// same physical quantities. If the user wants to shoot two particles with
// different momentum, position, etc., should invoke GeneratePrimaryVertex()
// method twice and set quantities on demand to the particle gun.

// Author: Makoto Asai, 1997
// --------------------------------------------------------------------
#ifndef G4ParticleGun_hh
#define G4ParticleGun_hh 1

#include "globals.hh"
#include "G4VPrimaryGenerator.hh"
#include "G4ThreeVector.hh"
#include "G4ParticleDefinition.hh"
#include "G4PrimaryVertex.hh"
#include "G4ParticleMomentum.hh"

class G4Event;
class G4ParticleGunMessenger;

class G4ParticleGun : public G4VPrimaryGenerator
{
  public:

    G4ParticleGun();
    explicit G4ParticleGun(G4int numberofparticles);
    explicit G4ParticleGun(G4ParticleDefinition* particleDef, 
                  G4int numberofparticles = 1);
      // Costructors. "numberofparticles" is the number of particles to be
      // shot at one invokation of GeneratePrimaryVertex() method.
      // All particles are shot with the same physical quantities.

    ~G4ParticleGun() override;

    G4ParticleGun(const G4ParticleGun&) = delete;
    const G4ParticleGun& operator=(const G4ParticleGun&) = delete;
    G4bool operator==(const G4ParticleGun&) const = delete;
    G4bool operator!=(const G4ParticleGun&) const = delete;

    void GeneratePrimaryVertex(G4Event* evt) override;
      // Creates a primary vertex at the given point
      // and put primary particles to it.

    // Followings are the Set methods for the particle properties.
    // SetParticleDefinition() should be called first.  
    // By using SetParticleMomentum(), both particle_momentum_direction and
    // particle_energy(Kinetic Energy) are set.
    //
    void SetParticleDefinition(G4ParticleDefinition* aParticleDefinition);
    void SetParticleEnergy(G4double aKineticEnergy);
    void SetParticleMomentum(G4double aMomentum);
    void SetParticleMomentum(G4ParticleMomentum aMomentum);
    inline void SetParticleMomentumDirection(G4ParticleMomentum aMomDirection)
      { particle_momentum_direction =  aMomDirection.unit(); }
    inline void SetParticleCharge(G4double aCharge)
      { particle_charge = aCharge; }
    inline void SetParticlePolarization(G4ThreeVector aVal)
      { particle_polarization = aVal; }
    inline void SetNumberOfParticles(G4int i)
      { NumberOfParticlesToBeGenerated = i; }

    inline G4ParticleDefinition* GetParticleDefinition() const
      { return particle_definition; }
    inline G4ParticleMomentum GetParticleMomentumDirection() const
      { return particle_momentum_direction; }
    inline G4double GetParticleEnergy() const
      { return particle_energy; }
    inline G4double GetParticleMomentum() const
      { return particle_momentum; }
    inline G4double GetParticleCharge() const
      { return particle_charge; }
    inline G4ThreeVector GetParticlePolarization() const
      { return particle_polarization; }
    inline G4int GetNumberOfParticles() const
      { return NumberOfParticlesToBeGenerated; }

  protected:  

     virtual void SetInitialValues();

     G4int                 NumberOfParticlesToBeGenerated = 0;
     G4ParticleDefinition* particle_definition = nullptr;
     G4ParticleMomentum    particle_momentum_direction;
     G4double              particle_energy = 0.0;
     G4double              particle_momentum = 0.0;
     G4double              particle_charge = 0.0;
     G4ThreeVector         particle_polarization;

  private:

     G4ParticleGunMessenger* theMessenger = nullptr;
};

#endif


