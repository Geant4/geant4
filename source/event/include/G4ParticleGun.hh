// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4ParticleGun.hh,v 1.4 2000-10-18 12:41:21 kurasige Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//

#ifndef G4ParticleGun_h
#define G4ParticleGun_h 1


#include "globals.hh"
#include "G4VPrimaryGenerator.hh"
#include "G4ThreeVector.hh"
#include "G4ParticleDefinition.hh"
#include "G4PrimaryVertex.hh"
#include "G4ParticleMomentum.hh"

class G4Event;
class G4ParticleGunMessenger;

// class description:
//
//  This is a concrete class of G4VPrimaryGenerator. It shoots a particle of given type
// to a given direction with a given kinetic energy. 
//  The position and time of the primary particle must be set by the corresponding
// set methods of G4VPrimaryGenerator base class, otherwise zero will be set.
//
//  The FAQ to this class is for randomizing position/direction/kinetic energy of primary
// particle. But, G4ParticleGun does NOT have any way of randomization. Instead, the user's
// concrete implementation of G4VUserPrimaryGeneratorAction which transmits G4Event object
// to this particle gun can randomize these quantities and set to this particle gun before
// invoking GeneratePrimaryVertex() method.
//  Note that, even if the particle gun shoots more than one particles at one invokation of
// GeneratePrimaryVertex() method, all particles have the same physical quantities. If the
// user wants to shoot two particles with different momentum, position, etc., invoke
// GeneratePrimaryVertex() method twice and set quantities on demand to the particle gun.
//

class G4ParticleGun:public G4VPrimaryGenerator
{
  public: // with description
     G4ParticleGun();
     G4ParticleGun(G4int numberofparticles);
     G4ParticleGun(G4ParticleDefinition * particleDef, 
                   G4int numberofparticles = 1);
     // costructors. "numberofparticles" is number of particles to be shoot at one invokation
     // of GeneratePrimaryVertex() method. All paricles are shoot with the same physical
     // quantities.

  public:
     virtual ~G4ParticleGun();
     G4ParticleGun(const G4ParticleGun &right);

     const G4ParticleGun & operator=(const G4ParticleGun &right);
     int operator==(const G4ParticleGun &right) const;
     int operator!=(const G4ParticleGun &right) const;

  public: // with description
     virtual void GeneratePrimaryVertex(G4Event* evt);
     // Creates a primary vertex at the given point and put primary particles to it.
     // Followings are set methods for the particle properties.
     //   SetParticleDefinition should be called first.  
     //   By using SetParticleMomentum(), both particle_momentum_direction and
     //   particle_energy(Kinetic Energy) are set.
     //   
     void SetParticleDefinition
       (G4ParticleDefinition * aParticleDefinition);
     void SetParticleMomentum(G4ParticleMomentum aMomentum);
     inline void SetParticleMomentumDirection
                 (G4ParticleMomentum aMomentumDirection)
     { particle_momentum_direction =  aMomentumDirection.unit(); }
     inline void SetParticleEnergy(G4double aKineticEnergy)
     { particle_energy = aKineticEnergy; }
     inline void SetParticlePosition(G4ThreeVector aPosition)
     { particle_position = aPosition; }
     inline void SetParticleTime(G4double aTime)
     { particle_time = aTime; }
     inline void SetParticleCharge(G4double aCharge)
     { particle_charge = aCharge; }
     inline void SetParticlePolarization(G4ThreeVector aVal)
     { particle_polarization = aVal; }
     inline void SetNumberOfParticles(G4int i)
     { NumberOfParticlesToBeGenerated = i; }

  public:
     inline G4ParticleDefinition* GetParticleDefinition()
     { return particle_definition; }
     inline G4ParticleMomentum GetParticleMomentumDirection()
     { return particle_momentum_direction; }
     inline G4double GetParticleEnergy()
     { return particle_energy; }
     inline G4double GetParticleCharge()
     { return particle_charge; }
     inline G4ThreeVector GetParticlePosition()
     { return particle_position; }
     inline G4double GetParticleTime()
     { return particle_time; }
     inline G4ThreeVector GetParticlePolarization()
     { return particle_polarization; }
     inline G4int GetNumberOfParticles()
     { return NumberOfParticlesToBeGenerated; }

  protected:  
     virtual void SetInitialValues();

     G4int                 NumberOfParticlesToBeGenerated;
     G4ParticleDefinition* particle_definition;
     G4ParticleMomentum    particle_momentum_direction;
     G4double	           particle_energy;
     G4double	           particle_charge;
     G4ThreeVector         particle_position;
     G4double              particle_time;
     G4ThreeVector         particle_polarization;

  private:
     G4ParticleGunMessenger* theMessenger;
};

#endif







