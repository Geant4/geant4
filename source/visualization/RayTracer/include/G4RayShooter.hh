// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4RayShooter.hh,v 1.2 2000-03-09 17:38:33 asaim Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//

// class description:
//
//  the object of this class shoots a ray (actually geantino) primary particle
// associated with a G4Event object. This slass must be used exclusively by
// G4RayTracer.
//

#ifndef G4RayShooter_h
#define G4RayShooter_h 1


#include "globals.hh"
#include "G4ThreeVector.hh"
#include "G4ParticleDefinition.hh"
#include "G4PrimaryVertex.hh"
#include "G4ParticleMomentum.hh"

class G4Event;

class G4RayShooter
{
  public: 
     G4RayShooter();
  public:
     virtual ~G4RayShooter();

  public: // with description
     void Shoot(G4Event* evt,G4ThreeVector vtx,G4ThreeVector direc);
     // This method generates a primary vertex and a primary particle at the
     // given vertex point and with the given direction. This method is invoked
     // by G4RayTracer.

  private:  
     void SetInitialValues();

     G4ParticleDefinition* particle_definition;
     G4ParticleMomentum    particle_momentum_direction;
     G4double	           particle_energy;
     G4ThreeVector         particle_position;
     G4double              particle_time;
     G4ThreeVector         particle_polarization;
};

#endif







