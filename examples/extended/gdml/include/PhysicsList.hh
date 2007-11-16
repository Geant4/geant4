#ifndef _PHYSICSLIST_H_
#define _PHYSICSLIST_H_

#include "G4VUserPhysicsList.hh"
#include "G4ParticleTypes.hh"
#include "globals.hh"

class PhysicsList: public G4VUserPhysicsList {
public:
   PhysicsList();
   ~PhysicsList();
protected:
   void ConstructParticle();
   void ConstructProcess();
   void SetCuts();
};

#endif







