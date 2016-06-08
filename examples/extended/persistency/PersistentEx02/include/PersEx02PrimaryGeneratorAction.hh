// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: PersEx02PrimaryGeneratorAction.hh,v 1.3 1999/11/29 18:23:32 morita Exp $
// GEANT4 tag $Name: geant4-01-01 $
//

#ifndef PersEx02PrimaryGeneratorAction_h
#define PersEx02PrimaryGeneratorAction_h 1

#include "G4VUserPrimaryGeneratorAction.hh"

class G4ParticleGun;
class G4Event;

class PersEx02PrimaryGeneratorAction : public G4VUserPrimaryGeneratorAction
{
  public:
    PersEx02PrimaryGeneratorAction();
    ~PersEx02PrimaryGeneratorAction();

  public:
    void GeneratePrimaries(G4Event* anEvent);

  private:
    G4ParticleGun* particleGun;
};

#endif


