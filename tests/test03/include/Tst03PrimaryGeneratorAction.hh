// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//

#ifndef Tst03PrimaryGeneratorAction_h
#define Tst03PrimaryGeneratorAction_h 1

#include "G4VUserPrimaryGeneratorAction.hh"

class G4ParticleGun;
class G4Event;

class Tst03PrimaryGeneratorAction : public G4VUserPrimaryGeneratorAction
{
  public:

    Tst03PrimaryGeneratorAction();
    ~Tst03PrimaryGeneratorAction();

  public:

    void GeneratePrimaries(G4Event* anEvent);

  private:

    G4ParticleGun* particleGun;
};

#endif /*Tst03PrimaryGeneratorAction_h*/
