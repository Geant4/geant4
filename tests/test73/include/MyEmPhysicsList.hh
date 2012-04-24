// User physics list
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#ifndef MyEmPhysicsList_h
#define MyEmPhysicsList_h 1

#include "G4VPhysicsConstructor.hh"
#include "globals.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class MyEmPhysicsList : public G4VPhysicsConstructor
{
	public: 
		MyEmPhysicsList(const G4String& name = "standard");
		~MyEmPhysicsList();

	public: 
		// This method is dummy for physics
		void ConstructParticle() {};
		// This method will be invoked in the Construct() method.
		// each physics process will be instantiated and
		// registered to the process manager of each particle type 
		void ConstructProcess();
		// Sets the step limit based on the logical volume
		virtual void AddStepMax(G4ParticleDefinition* particle,	G4ProcessManager* pmanager); 

	private:
		G4int  verbose;

};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif








