#ifndef BrachyFactory_h
#define BrachyFactory_h 1
#include "G4VUserPrimaryGeneratorAction.hh"
#include "G4RadioactiveDecay.hh"
class G4ParticleGun;
class G4Run;
class G4Event;
class BrachyAnalysisManager;

class BrachyFactory
{
public:
  BrachyFactory();
 virtual ~BrachyFactory();
  virtual  G4VUserPrimaryGeneratorAction* CreatePrimaryGeneratorAction()=0;
  virtual void CreateSource(G4VPhysicalVolume*)=0;
  virtual void CleanSource()=0;
};
#endif

















