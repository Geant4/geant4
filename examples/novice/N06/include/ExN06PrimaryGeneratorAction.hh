// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//

#ifndef ExN06PrimaryGeneratorAction_h
#define ExN06PrimaryGeneratorAction_h 1

#include "G4VUserPrimaryGeneratorAction.hh"

class G4ParticleGun;
class G4Event;

class ExN06PrimaryGeneratorAction : public G4VUserPrimaryGeneratorAction
{
  public:

    ExN06PrimaryGeneratorAction();
    ~ExN06PrimaryGeneratorAction();

  public:

    void GeneratePrimaries(G4Event* anEvent);

  private:

    G4ParticleGun* particleGun;
};

#endif /*ExN06PrimaryGeneratorAction_h*/
