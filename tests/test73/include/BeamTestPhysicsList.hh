// modular list on base of old PhysicsList
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#ifndef PhysicsList_h
#define PhysicsList_h 1

#include "G4VModularPhysicsList.hh"
#include "globals.hh"

class G4VPhysicsConstructor;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class BeamTestPhysicsList: public G4VModularPhysicsList
{
	public:
		BeamTestPhysicsList();
		~BeamTestPhysicsList();

		void ConstructParticle();
		void ConstructProcess();
		void AddPhysicsList();
 		// Sets the step limit based on the logical volume
		//void AddStepMax(); 
		void SetCuts();
		//void SetCutForGamma(G4double);
		//void SetCutForElectron(G4double);
		//void SetCutForPositron(G4double);

	private:
		G4double cutForGamma;
		G4double cutForElectron;
		G4double cutForPositron;
		G4double cutForMuon;
		G4double cutForPion;
		//G4double cutForPositron;
		G4double currentDefaultCut;

		G4VPhysicsConstructor*  emPhysicsList;
		G4String name;

};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif

