// Created by Joanna Weng 9.11.2004

#include "G4Electron.hh"
#include "G4Positron.hh"

#include "GFlashParticleBounds.hh"

GFlashParticleBounds::GFlashParticleBounds()
{    
	// e+e- defaults
	EMinEneToParametrise = 0.10*GeV;  
	EMaxEneToParametrise = 10000.00*GeV;  
	EEneToKill = 0.1*GeV; // Energie at which electrons are killed
}

GFlashParticleBounds::~GFlashParticleBounds()
{
}

void GFlashParticleBounds::SetMinEneToParametrise(G4ParticleDefinition &particleType,G4double enemin)
{ 
	if( &particleType == G4Electron::ElectronDefinition()||
		&particleType == G4Positron::PositronDefinition()) 
	EMinEneToParametrise = enemin;
}

void GFlashParticleBounds::SetMaxEneToParametrise(G4ParticleDefinition &particleType,G4double enemax)
{
	if( &particleType == G4Electron::ElectronDefinition()||
		&particleType == G4Positron::PositronDefinition()) 
	EMaxEneToParametrise = enemax; 
}

void GFlashParticleBounds::SetEneToKill(G4ParticleDefinition &particleType,G4double enekill)
{
	if( &particleType == G4Electron::ElectronDefinition()||
		&particleType == G4Positron::PositronDefinition()) 
	EEneToKill = enekill; 
}

G4double GFlashParticleBounds::GetMinEneToParametrise(G4ParticleDefinition &particleType) 
{ 
	G4double result = DBL_MAX;
	if( &particleType == G4Electron::ElectronDefinition()||
		&particleType == G4Positron::PositronDefinition()) 
	{
		result = EMinEneToParametrise;
	}        
	return result;
}

G4double GFlashParticleBounds::GetMaxEneToParametrise(G4ParticleDefinition &particleType) 
{ 
	G4double result = 0;
	if( &particleType == G4Electron::ElectronDefinition()||
		&particleType == G4Positron::PositronDefinition()) 
	{
		result = EMaxEneToParametrise; 
	}
	return result;
}

G4double GFlashParticleBounds::GetEneToKill(G4ParticleDefinition & particleType) 
{
	if (&particleType == G4Electron::ElectronDefinition() ||
		&particleType == G4Positron::PositronDefinition())  
	return EEneToKill;
	else return (-DBL_MAX);
}





