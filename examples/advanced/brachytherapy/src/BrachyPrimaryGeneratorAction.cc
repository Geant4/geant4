//    ********************************************
//    *                                          *
//    *    BrachyPrimaryGeneratorAction.cc       *
//    *                                          *
//    ********************************************

#include "BrachyPrimaryGeneratorAction.hh"

#include "G4ParticleTable.hh"
#include "Randomize.hh"  
#include "G4Event.hh"
#include "G4ParticleGun.hh"
#include "G4IonTable.hh"
#include "G4RadioactiveDecay.hh"
#include "G4UImanager.hh"
#include "globals.hh"
#include <math.h>

//....

BrachyPrimaryGeneratorAction::BrachyPrimaryGeneratorAction()
{
 // Generate a gamma particle with energy = Ir-192 mean energy
 G4int NumParticles = 1;
 G4double Energy = 0.356*MeV;

 m_pParticleGun = new G4ParticleGun(NumParticles);
 if(m_pParticleGun) 
	m_pParticleGun->SetParticleEnergy(Energy);
}

//....

BrachyPrimaryGeneratorAction::~BrachyPrimaryGeneratorAction()
{
 if(m_pParticleGun)
	delete m_pParticleGun;
}

//....

void BrachyPrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent)
{
 G4ParticleTable* pParticleTable = G4ParticleTable::GetParticleTable();
 G4String ParticleName = "gamma";
 G4ParticleDefinition* pParticle = pParticleTable->FindParticle(ParticleName);

 m_pParticleGun->SetParticleDefinition(pParticle);

 // Random generation of gamma source point inside the Iridium core cylinder(R=0.3*mm,h=3.5*mm)
 G4double radius = 0.3*mm;
 G4double x,y,z;
 do{
   x = (G4UniformRand()-0.5)*radius/0.5;
   y = (G4UniformRand()-0.5)*radius/0.5;
   }while(x*x+y*y > radius*radius);
 z = (G4UniformRand()-0.5)*1.75*mm/0.5;    

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

 m_pParticleGun->GeneratePrimaryVertex(anEvent);
}
