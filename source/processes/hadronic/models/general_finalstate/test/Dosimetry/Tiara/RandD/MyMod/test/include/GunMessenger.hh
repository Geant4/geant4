#ifndef __GUNMESSENGER_DEFINED__
#define __GUNMESSENGER_DEFINED__

#include "G4UImessenger.hh"
#include "globals.hh"

class G4UIcmdWithAString;
class G4UIcmdWithoutParameter;
class ParticleGun;
class G4UIcommand;

class GunMessenger : public G4UImessenger
{
public:
  GunMessenger(ParticleGun* pGun);
  ~GunMessenger();
  void SetNewValue(G4UIcommand* pCmd,G4String szValue);
private:
  G4UIcmdWithAString*      m_pCmdSetNeutronGun;
  G4UIcmdWithoutParameter* m_pCmdSetProtonGun;

  ParticleGun*             m_pGun;
};


#endif
