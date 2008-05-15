
#include "PrimaryGeneratorAction.hh"

//Geant4 includes
#include "globals.hh"
#include "G4ParticleDefinition.hh"

PrimaryGeneratorAction::PrimaryGeneratorAction()
{
	//particle gun and particle table 
	particleGun = new G4ParticleGun();
  	particleTable = G4ParticleTable::GetParticleTable();

	//default particle
  	particleGun->SetParticleDefinition(particleTable->FindParticle("geantino"));
  	particleGun->SetParticleEnergy( 1.0*MeV );
 
	// particleGun->SetParticlePosition(G4ThreeVector(-140,1100,-160));
        // error from SBT test for Polycone
        G4ThreeVector err1=G4ThreeVector(-1260,-560,40);//outside
        G4ThreeVector err2=G4ThreeVector(100,-240,120);//inside
        G4ThreeVector err2v=(err2-err1).unit();
  
 particleGun->SetParticleMomentumDirection(err2v);

 
  particleGun->SetParticlePosition(err1);


} //end of constructor

//destructor
PrimaryGeneratorAction::~PrimaryGeneratorAction()
{
  	delete particleGun;
}


void PrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent)
{
  	particleGun->GeneratePrimaryVertex(anEvent);
}

//EOF PrimaryGeneratorAction.cc
