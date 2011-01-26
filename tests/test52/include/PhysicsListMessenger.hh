#ifndef PHYSICSLISTMESSENGER_HH
#define PHYSICSLISTMESSENGER_HH 1

#include "G4UImessenger.hh"
#include "globals.hh"

class PhysicsList;
class G4UIdirectory;
class G4UIcmdWithAString;
class G4UIcmdWithADoubleAndUnit;

class PhysicsListMessenger : public G4UImessenger {

 public:
   PhysicsListMessenger(PhysicsList*);
   ~PhysicsListMessenger();

   void SetNewValue(G4UIcommand*, G4String);

 private:
   PhysicsList* physicsList;

   G4UIdirectory* physicsDirectory;
   G4UIcmdWithAString* physicsConstrCmd;
   G4UIcmdWithADoubleAndUnit* prodThresholdCmd;
};

#endif // PHYSICSLISTMESSENGER_HH
