//
// ********************************************************************
// * DISCLAIMER                                                       *
// *                                                                  *
// * The following disclaimer summarizes all the specific disclaimers *
// * of contributors to this software. The specific disclaimers,which *
// * govern, are listed with their locations in:                      *
// *   http://cern.ch/geant4/license                                  *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.                                                             *
// *                                                                  *
// * This  code  implementation is the  intellectual property  of the *
// * GEANT4 collaboration.                                            *
// * By copying,  distributing  or modifying the Program (or any work *
// * based  on  the Program)  you indicate  your  acceptance of  this *
// * statement, and all its terms.                                    *
// ********************************************************************
//
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
