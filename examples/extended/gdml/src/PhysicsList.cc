#include "PhysicsList.hh"

PhysicsList::PhysicsList() {
}

PhysicsList::~PhysicsList() {
}

void PhysicsList::ConstructParticle() {

   G4Geantino::GeantinoDefinition();
}

void PhysicsList::ConstructProcess() {

   AddTransportation();
}

void PhysicsList::SetCuts() {

   G4int temp = GetVerboseLevel();
   SetVerboseLevel(0);                                                           
   SetCutsWithDefault();   
   SetVerboseLevel(temp);  
}

