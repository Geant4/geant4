#ifndef __PARTICLE_GUN_DEFINED__
#define __PARTICLE_GUN_DEFINED__

#include "G4VUserPrimaryGeneratorAction.hh"
#include "globals.hh"
#include "G4ParticleGun.hh"

class G4Event;
class IVector;
class GunMessenger;

class ParticleGun : public G4VUserPrimaryGeneratorAction
{
public:
  ParticleGun();
  ~ParticleGun();
  void GeneratePrimaries(G4Event* pEvent);
  G4double CurrEnergy(){return pGun->GetParticleEnergy();};
  void CalculateNeutrons(bool bCalculate,char* szFileName);
  bool GetStat(){return m_bCalculateNeutrons;};
  G4double GetMax();
  G4double GetArea(){return m_dArea;}
  G4double GetRArea(){return m_dRArea;}
private:
  G4double GetEnergy();

  bool           m_bCalculateNeutrons;
  G4double       m_dArea;
  G4double       m_dRArea;
  IVector*       m_pNeutronData;
  G4ParticleGun* pGun;
  GunMessenger*  m_pMessenger;
};
#endif
