//
// FredPhysicsList
//
// Implementation of Physics for fred
//

#include "FredPhysicsList.hh"
#include "G4ParticleTypes.hh"

//
// Constructor
//
FredPhysicsList::FredPhysicsList(): G4VUserPhysicsList()
{
  defaultCutValue = 2.0*mm;
}



FredPhysicsList::~FredPhysicsList() {;}

//
// ConstructParticle: Define particles
//
void FredPhysicsList::ConstructParticle()
{
	// Particles are defined via static member functions

	G4Geantino::GeantinoDefinition();
	G4MuonMinus::MuonMinusDefinition();
}

//
// ConstructProcess: define processes
//
void FredPhysicsList::ConstructProcess()
{
	AddTransportation();
}

//
// SetCuts: define cuts
//
void FredPhysicsList::SetCuts()
{
	// From "novice" example 1

	// uppress error messages even in case e/gamma/proton do not exist            
	G4int temp = GetVerboseLevel();                                                SetVerboseLevel(0);                                                           
	//  " G4VUserPhysicsList::SetCutsWithDefault" method sets 
	//   the default cut value for all particle types 
	SetCutsWithDefault();   

	// Retrieve verbose level
	SetVerboseLevel(temp);  
	//  " G4VUserPhysicsList::SetCutsWithDefault" method sets 
	//   the default cut value for all particle types 
	SetCutsWithDefault();   
}
