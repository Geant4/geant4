#include "BrachyFactoryIr.hh"
#include"BrachyPrimaryGeneratorActionIr.hh"
#include "G4ParticleTable.hh"
#include "Randomize.hh"  
#include "G4Event.hh"
#include "G4ParticleGun.hh"
#include "G4IonTable.hh"
#include "G4RadioactiveDecay.hh"
#include "G4UImanager.hh"
#include "globals.hh"
#include <math.h>
#include "G4RunManager.hh" 
#include"BrachyDetectorMessenger.hh"
#include"BrachyDetectorConstructionIr.hh"
BrachyFactoryIr:: BrachyFactoryIr()
{
   pIridio=new  BrachyDetectorConstructionIr();

}

BrachyFactoryIr:: ~BrachyFactoryIr()
{
  delete pIridio;

}
 

G4VUserPrimaryGeneratorAction*  BrachyFactoryIr::CreatePrimaryGeneratorAction()

{ 
 
   
  G4VUserPrimaryGeneratorAction*    pIridium =new BrachyPrimaryGeneratorActionIr();
 if(pIridium) return pIridium ;


}

void BrachyFactoryIr::CreateSource(G4VPhysicalVolume* mother)
{

  pIridio -> ConstructIridium(mother);

}
void BrachyFactoryIr::CleanSource()
{

  pIridio -> CleanIridium();

}
