//    ********************************************
//    *                                          *
//    *    BrachyPrimaryGeneratorAction.cc       *
//    *                                          *
//    ********************************************

#include "BrachyPrimaryGeneratorActionI.hh"
#include "BrachyAnalysisManager.hh"
//#include "BrachyPrimaryGeneratorMessenger.hh"
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

//....

BrachyPrimaryGeneratorActionI::BrachyPrimaryGeneratorActionI()
{
 
 G4int NumParticles = 1;
 
  m_pParticleGun = new G4ParticleGun(NumParticles);

//create a messenger for this class
  // gunMessenger = new BrachyPrimaryGeneratorMessenger(this);

   vettore.push_back(0.783913);
  vettore.push_back(0.170416);
   vettore.push_back(0.045671);
 
}

//....

BrachyPrimaryGeneratorActionI::~BrachyPrimaryGeneratorActionI()
{
 if(m_pParticleGun)
	delete m_pParticleGun;

 //delete gunMessenger;
}

//....

void BrachyPrimaryGeneratorActionI::GeneratePrimaries(G4Event* anEvent)
{
BrachyAnalysisManager* analysis = BrachyAnalysisManager::getInstance();

 G4ParticleTable* pParticleTable = G4ParticleTable::GetParticleTable();
 G4String ParticleName = "gamma";
 G4ParticleDefinition* pParticle = pParticleTable->FindParticle(ParticleName);
 // G4cout<<"arrivo allfind particle"<<G4endl;

 m_pParticleGun->SetParticleDefinition(pParticle);

 //  Random generation of gamma source point inside the Iodium core
 G4double x,y,z;
 G4double radiusmax = 0.30*mm;
 G4double radiusmin=0.085*mm;
 
 do{
   x = (G4UniformRand()-0.5)*(radiusmax)/0.5;
   y = (G4UniformRand()-0.5)*(radiusmax)/0.5;
   }while(((x*x+y*y )> (radiusmax*radiusmax))||((x*x+y*y)<(radiusmin*radiusmin)));
 
 z = (G4UniformRand()-0.5)*1.75*mm/0.5 ;

 G4ThreeVector position(x,y,z);
 m_pParticleGun->SetParticlePosition(position);


 // Random generation of the impulse direction
 G4double a,b,c;
 G4double n;
 do{
   a = (G4UniformRand()-0.5)/0.5;
   b = (G4UniformRand()-0.5)/0.5; 
   c = (G4UniformRand()-0.5)/0.5;
   n = a*a+b*b+c*c;
   }while(n > 1 || n == 0.0);
 n = sqrt(n);
 a /= n;
 b /= n;
 c /= n;

 G4ThreeVector direction(a,b,c);
 m_pParticleGun->SetParticleMomentumDirection(direction);


 
 G4double random=G4UniformRand();
 G4double sum=0;
 G4int i=0;
 while(sum<random){sum+=vettore[i];
                       i++;}

 

 if(i==1){Energy=27.4*keV;}
 else{ 
       if(i==2){Energy=31.4*keV;}
       else {Energy=35.5*keV;}}
	       
      //G4double Energy= 28.451602*keV;

   m_pParticleGun->SetParticleEnergy(Energy);
   
   analysis->Spectrum(Energy);
  
   m_pParticleGun->GeneratePrimaryVertex(anEvent);
 
}
 

G4double  BrachyPrimaryGeneratorActionI::GetEnergy()
{ return Energy;}

