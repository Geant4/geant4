#include "BrachyFactoryI.hh"
#include"BrachyPrimaryGeneratorActionI.hh"
#include"BrachyDetectorConstructionI.hh"
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

BrachyFactoryI:: BrachyFactoryI()
{
   pIodio=new  BrachyDetectorConstructionI();
}

BrachyFactoryI:: ~BrachyFactoryI()
{
 
delete pIodio;
}
 

G4VUserPrimaryGeneratorAction*  BrachyFactoryI::CreatePrimaryGeneratorAction()

{ 
  
   
  G4VUserPrimaryGeneratorAction*    pIodium =new BrachyPrimaryGeneratorActionI();
 if(pIodium) return pIodium ;


 

}

void BrachyFactoryI::CreateSource(G4VPhysicalVolume* mother)
{

  pIodio -> ConstructIodium(mother);

}

void BrachyFactoryI::CleanSource()
{;}

