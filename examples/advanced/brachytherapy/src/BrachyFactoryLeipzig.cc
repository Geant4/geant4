#include "BrachyFactoryLeipzig.hh"
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
#include"BrachyDetectorConstructionLeipzig.hh"
BrachyFactoryLeipzig:: BrachyFactoryLeipzig()
{
   pLeipzig=new  BrachyDetectorConstructionLeipzig();
 
}

BrachyFactoryLeipzig:: ~BrachyFactoryLeipzig()
{
  delete pLeipzig;
 
}
 

G4VUserPrimaryGeneratorAction*  BrachyFactoryLeipzig::CreatePrimaryGeneratorAction()

{ 
 
   
  G4VUserPrimaryGeneratorAction*    pIridium =new BrachyPrimaryGeneratorActionIr();
 if(pIridium) return pIridium ;


}

void BrachyFactoryLeipzig::CreateSource(G4VPhysicalVolume* mother)
{

  pLeipzig -> ConstructLeipzig(mother);

}
void BrachyFactoryLeipzig::CleanSource()
{

  ;

}
