
#ifndef PrimaryGeneratorAction_h
#define PrimaryGeneratorAction_h 1

//Geant4 classes
#include "G4ParticleTable.hh"
#include "G4ParticleGun.hh"
#include "G4Event.hh"

//Virtual class
#include "G4VUserPrimaryGeneratorAction.hh"

class PrimaryGeneratorAction : public G4VUserPrimaryGeneratorAction
{
  public:
	//constructor and destructor
    	PrimaryGeneratorAction();
    	~PrimaryGeneratorAction();

	//used by Geant4 to generate the primary particles of the event
    	void GeneratePrimaries(G4Event* anEvent);

  private:

    	G4ParticleGun* particleGun;
    	G4ParticleTable * particleTable;

};

#endif

//EOF PrimaryGeneratorAction.hh
