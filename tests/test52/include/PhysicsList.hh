#ifndef PHYSICSLIST_HH
#define PHYSICSLIST_HH 1

#include "G4VModularPhysicsList.hh"
#include "globals.hh"

class PhysicsListMessenger;


class PhysicsList : public G4VModularPhysicsList {

 public:
   PhysicsList();
   virtual ~PhysicsList();

   void RegisterPhysConstructor(const G4String& constrName);
   void SetProdThreshold(double cut);
   void SetCuts();

 private:
   PhysicsListMessenger* messenger;

   static G4VPhysicsConstructor* emElectron;
   static G4VPhysicsConstructor* emPositron;
   static G4VPhysicsConstructor* emPhoton;
};

#endif // PHYSICSLIST_HH
