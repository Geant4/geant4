// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4ParticleGun.hh,v 1.1 1999-01-07 16:06:33 gunter Exp $
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

class G4ParticleGun:public G4VPrimaryGenerator
{
  public:
     G4ParticleGun();
     G4ParticleGun(G4int numberofparticles);
     G4ParticleGun(G4ParticleDefinition * particleDef, 
                   G4int numberofparticles = 1);
     virtual ~G4ParticleGun();
     G4ParticleGun(const G4ParticleGun &right);

     const G4ParticleGun & operator=(const G4ParticleGun &right);
     int operator==(const G4ParticleGun &right) const;
     int operator!=(const G4ParticleGun &right) const;

     virtual void GeneratePrimaryVertex(G4Event* evt);
  
     // Set particle properties to be shoot out
     //   SetParticleDefinition should be called at first  
     //   By using SetParticleMomentum(), both particle_momentum_direction and
     //   particle_energy(Kinetic Energy) are determined.

     void SetParticleMomentum(G4ParticleMomentum aMomentum);
     inline void SetParticleDefinition
                 (G4ParticleDefinition * aParticleDefinition)
     { particle_definition = aParticleDefinition; }
     inline void SetParticleMomentumDirection
                 (G4ParticleMomentum aMomentumDirection)
     { particle_momentum_direction =  aMomentumDirection.unit(); }
     inline void SetParticleEnergy(G4double aKineticEnergy)
     { particle_energy = aKineticEnergy; }
     inline void SetParticlePosition(G4ThreeVector aPosition)
     { particle_position = aPosition; }
     inline void SetParticleTime(G4double aTime)
     { particle_time = aTime; }
     inline G4ParticleDefinition* GetParticleDefinition()
     { return particle_definition; }
     inline G4ParticleMomentum GetParticleMomentumDirection()
     { return particle_momentum_direction; }
     inline G4double GetParticleEnergy()
     { return particle_energy; }
     inline G4ThreeVector GetParticlePosition()
     { return particle_position; }
     inline G4double GetParticleTime()
     { return particle_time; }
     inline void SetParticlePolarization(G4ThreeVector aVal)
     { particle_polarization = aVal; }
     inline G4ThreeVector GetParticlePolarization()
     { return particle_polarization; }
     inline void SetNumberOfParticles(G4int i)
     { NumberOfParticlesToBeGenerated = i; }
     inline G4int GetNumberOfParticles()
     { return NumberOfParticlesToBeGenerated; }

  protected:  
     virtual void SetInitialValues();

     G4int                 NumberOfParticlesToBeGenerated;
     G4ParticleDefinition* particle_definition;
     G4ParticleMomentum    particle_momentum_direction;
     G4double	           particle_energy;
     G4ThreeVector         particle_position;
     G4double              particle_time;
     G4ThreeVector         particle_polarization;

  private:
     G4ParticleGunMessenger* theMessenger;
};

#endif







