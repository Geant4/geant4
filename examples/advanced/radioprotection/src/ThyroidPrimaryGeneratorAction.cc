//
//
//    ********************************************
//    *                                          *
//    *    ThyroidPrimaryGeneratorAction.cc      *
//    *                                          *
//    ********************************************

#include "ThyroidPrimaryGeneratorAction.hh"

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

ThyroidPrimaryGeneratorAction::ThyroidPrimaryGeneratorAction()
{
 // Generate a gamma particle with energy = I-131 mean energy

    G4int NumParticles = 1; 
    m_pParticleGun = new G4ParticleGun(NumParticles);

    // G4double Energy = 0.364*MeV;

if(m_pParticleGun) 
       m_pParticleGun->SetParticleEnergy(Energy);
}

//....

ThyroidPrimaryGeneratorAction::~ThyroidPrimaryGeneratorAction()
  {
   if(m_pParticleGun)
  delete m_pParticleGun;
  }

//....

void ThyroidPrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent)
{
 G4ParticleTable* pParticleTable = G4ParticleTable::GetParticleTable();
 G4String ParticleName = "gamma";
 G4ParticleDefinition* pParticle = pParticleTable->FindParticle(ParticleName);

 m_pParticleGun->SetParticleDefinition(pParticle);

 // Random generation of gamma source point inside the Iodine Nodule
 G4double radius = 0.6*cm;
 G4double x,y,z;
 do{
   x = (G4UniformRand()-0.5)*radius/0.5;
   y = (G4UniformRand()-0.5)*radius/0.5;
   z = (G4UniformRand()-0.5)*radius/0.5; 
   }
   while(x*x+y*y+z*z > radius*radius);
 
  x = x-2.0;
 y = y-4.0;
  z = z-1.0;
     //
 G4ThreeVector position(x,y,z);
 m_pParticleGun->SetParticlePosition(position);

 // Random generation of the impulse direction
 G4double a,b,c;
 G4double n;
 G4double Energy;
 G4double Pr;
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
 Pr=G4UniformRand()*100;
 //
 // Gamma Spectrum in 131-I decay taken from http://www2.bnl.gov/ton  tables
 // Selected most relevant energies and renormalized to 1
 //

 if(Pr<2.62)                {Energy=0.080*MeV;}
 if(Pr>2.62  && Pr< 8.75)   {Energy=0.284*MeV;}
 if(Pr>8.75  && Pr<90.95)   {Energy=0.364*MeV;}
 if(Pr>90.95 && Pr<98.19)   {Energy=0.637*MeV;}
 if(Pr>98.19)               {Energy=0.723*MeV;}
 //
 m_pParticleGun->SetParticleMomentumDirection(direction); 
 m_pParticleGun->SetParticleEnergy(Energy);
 m_pParticleGun->GeneratePrimaryVertex(anEvent);

}







