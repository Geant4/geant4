// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4RayShooter.hh,v 1.1 2000-01-29 00:44:49 asaim Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
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

  public:
     void Shoot(G4Event* evt,G4ThreeVector vtx,G4ThreeVector direc);

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







