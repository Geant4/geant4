// Created by Joanna Weng, 9.11.04

#ifndef GFlashShowerModelMessenger_h
#define GFlashShowerModelMessenger_h 1

#include "globals.hh"
#include "G4UImessenger.hh"

class GFlashShowerModel;
class G4UIdirectory;
class G4UIcmdWithAString;
class G4UIcmdWithAnInteger;
class G4UIcmdWithADoubleAndUnit;
class G4UIcmdWithADouble;

class GFlashShowerModelMessenger: public G4UImessenger
{
	public:
	GFlashShowerModelMessenger(GFlashShowerModel * myModel);
	~GFlashShowerModelMessenger();
	
	void SetNewValue(G4UIcommand * command,G4String newValues);
	G4String GetCurrentValue(G4UIcommand * command);
	
	private:
	GFlashShowerModel* myModel;
	G4UIdirectory* 	myParaDir;
	G4UIcmdWithAString*	SwitchCmd;
	G4UIcmdWithAnInteger*	FlagCmd;
	G4UIcmdWithAnInteger*	ContCmd; // Containment Check
	G4UIcmdWithADouble* 	StepInX0Cmd;		
 	
	G4UIcmdWithADoubleAndUnit* 	EmaxCmd;
	G4UIcmdWithADoubleAndUnit* 	EminCmd;
	G4UIcmdWithADoubleAndUnit* 	EkillCmd;
};

#endif













