///////////////////////////////////////////////////////////////////////////////
// File: CMSPrimaryGeneratorMessenger.hh
// Author: I. Gonzalez (based on Geant4 examples)
// Description: Adds a new command to (un)select random shooting.
// Modification: 18/04/00  P.Arce New commands
///////////////////////////////////////////////////////////////////////////////

#ifndef CMSPrimaryGeneratorMessenger_h
#define CMSPrimaryGeneratorMessenger_h 1

#include "globals.hh"
#include "G4UImessenger.hh"

class CMSPrimaryGeneratorAction;

#include "G4UIcmdWithAString.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"
#include "G4UIcmdWithADouble.hh"
#include "G4UIcmdWithAnInteger.hh"

class CMSPrimaryGeneratorMessenger: public G4UImessenger {
public:
  CMSPrimaryGeneratorMessenger(CMSPrimaryGeneratorAction* myGun);
  ~CMSPrimaryGeneratorMessenger();
  
  void SetNewValue(G4UIcommand * command,G4String newValues);
  
private:
  CMSPrimaryGeneratorAction* myAction;
  
  G4UIcmdWithAnInteger*      verboseCmd;
  G4UIcmdWithAString*        rndmCmd;
  G4UIcmdWithAString*        scanCmd;
  G4UIcmdWithADoubleAndUnit* minEnergyCmd;
  G4UIcmdWithADoubleAndUnit* maxEnergyCmd;
  G4UIcmdWithADoubleAndUnit* minPhiCmd;
  G4UIcmdWithADoubleAndUnit* maxPhiCmd;
  G4UIcmdWithADouble*        minEtaCmd; 
  G4UIcmdWithADouble*        maxEtaCmd; 
  G4UIcmdWithAnInteger*      stepsPhiCmd;
  G4UIcmdWithAnInteger*      stepsEtaCmd;
  G4UIcmdWithAnInteger*      runNoCmd;
};

#endif

