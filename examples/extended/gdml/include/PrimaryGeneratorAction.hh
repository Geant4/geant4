#ifndef _PRIMARYGENERATORACTION_H_
#define _PRIMARYGENERATORACTION_H_

#include "G4VUserPrimaryGeneratorAction.hh"

#include "G4Event.hh"
#include "G4ParticleGun.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include "globals.hh"

//class G4ParticleGun;
//class G4Event;

class PrimaryGeneratorAction : public G4VUserPrimaryGeneratorAction {
public:
   PrimaryGeneratorAction();
   ~PrimaryGeneratorAction();

   void GeneratePrimaries(G4Event* anEvent);
private:
    G4ParticleGun* particleGun;
};

#endif


