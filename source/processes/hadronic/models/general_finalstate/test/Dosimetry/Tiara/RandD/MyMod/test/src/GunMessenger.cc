#include "GunMessenger.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIcmdWithoutParameter.hh"
#include "ParticleGun.hh"

//#include "stdafx.h"
GunMessenger::GunMessenger(ParticleGun* pGun)
{

  m_pGun = pGun;

  m_pCmdSetNeutronGun = new G4UIcmdWithAString("/gun/setNeutron",this);
  m_pCmdSetNeutronGun->SetGuidance("Sets a neutron event generator");
  m_pCmdSetNeutronGun->SetParameterName("DataFile",false);
  m_pCmdSetNeutronGun->AvailableForStates(PreInit,Idle);

  m_pCmdSetProtonGun = new G4UIcmdWithoutParameter("/gun/setProton",this);
  m_pCmdSetProtonGun->SetGuidance("Sets a proton(standard) event generator");
  m_pCmdSetProtonGun->AvailableForStates(PreInit,Idle);
}

GunMessenger::~GunMessenger()
{
  delete m_pCmdSetNeutronGun;
  delete m_pCmdSetProtonGun;
}

void GunMessenger::SetNewValue(G4UIcommand* pCmd,G4String szValue)
{
  if(pCmd == m_pCmdSetNeutronGun)
    m_pGun->CalculateNeutrons(true,(char*)szValue.data());
  else if(pCmd == m_pCmdSetProtonGun)
    m_pGun->CalculateNeutrons(false,NULL);
}
