#ifndef __PHYSICS_MESSENGER_DEFINED__
#define __PHYSICS_MESSENGER_DEFINED__

#include "G4UImessenger.hh"
#include "globals.hh"

class G4UIcmdWithAString;
class G4UIcommand;
class PhysicsList;

class G4DosimPhysicsMessenger : public G4UImessenger
{
public:
  G4DosimPhysicsMessenger(PhysicsList* pMother);
  ~G4DosimPhysicsMessenger();
  void SetNewValue(G4UIcommand* pCmd,G4String szValue);
private:
  G4UIcmdWithAString* m_pChangeModel;
  PhysicsList* m_pMother;
};
#endif
