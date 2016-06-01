// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4RayPrimaryGeneratorAction.hh,v 2.1 1998/07/12 03:09:19 urbi Exp $
// GEANT4 tag $Name: geant4-00 $
//

#ifndef G4RayPrimaryGeneratorAction_h
#define G4RayPrimaryGeneratorAction_h 1

#include "G4VUserPrimaryGeneratorAction.hh"

class G4ParticleGun;
class G4Event;

class G4RayPrimaryGeneratorAction : public G4VUserPrimaryGeneratorAction
{
  public:
    G4RayPrimaryGeneratorAction();
    ~G4RayPrimaryGeneratorAction();

  public:
    void GeneratePrimaries(G4Event* anEvent);
    G4ParticleGun* GetParticleGun();

  private:
    G4ParticleGun* particleGun;
};

#endif


