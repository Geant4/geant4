#ifndef BrachyFactoryI_h
#define BrachyFactoryI_h 1
#include "G4VUserPrimaryGeneratorAction.hh"
#include "G4RadioactiveDecay.hh"
#include"BrachyDetectorConstructionI.hh"
#include"BrachyFactory.hh"
#include "G4RunManager.hh"
class G4ParticleGun;
class G4Run;
class G4Event;
class BrachyAnalysisManager;
class BrachyFactory;
class BrachyPrimaryGeneratorActionI;
class BrachyDetectorConstructionI;
class BrachyFactoryI:public BrachyFactory
{
public:
  BrachyFactoryI();
 ~BrachyFactoryI();
  G4VUserPrimaryGeneratorAction* CreatePrimaryGeneratorAction();
  void CreateSource(G4VPhysicalVolume*);
 void CleanSource();
private:
  BrachyDetectorConstructionI* pIodio;
 

};
#endif
