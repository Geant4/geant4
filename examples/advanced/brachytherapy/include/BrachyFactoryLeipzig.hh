#ifndef BrachyFactoryLeipzig_h
#define BrachyFactoryLeipzig_h 1
#include "G4VUserPrimaryGeneratorAction.hh"
#include "G4RadioactiveDecay.hh"
#include"BrachyFactory.hh"
#include "G4RunManager.hh"
class G4ParticleGun;
class G4Run;
class G4Event;
class BrachyAnalysisManager;
class BrachyFactory;
class BrachyPrimaryGeneratorActionIr;
class BrachyDetectorConstructionLeipzig;

class BrachyFactoryLeipzig:public BrachyFactory
{
public:
  BrachyFactoryLeipzig();
 ~BrachyFactoryLeipzig();
  G4VUserPrimaryGeneratorAction* CreatePrimaryGeneratorAction();
  void CreateSource(G4VPhysicalVolume*);
  void CleanSource();
private:
  BrachyDetectorConstructionLeipzig* pLeipzig;
};
#endif
