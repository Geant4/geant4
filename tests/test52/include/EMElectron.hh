#ifndef EMELECTRON_HH
#define EMELECTRON_HH 1

#include "globals.hh"

class G4VPhysicsConstructor;


class EMElectron {

 public:
   static EMElectron* Instance();
   ~EMElectron();

   G4VPhysicsConstructor* PhysicsConstructor(const G4String& name="");

 protected: 
   EMElectron();

 private: 
   static EMElectron* instance;
   static G4VPhysicsConstructor* physConstructor;
};

#endif EMELECTRON_HH
