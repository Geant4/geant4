///////////////////////////////////////////////////////////////////////////////
// File: CCalPrimaryGeneratorMessenger.hh
// Description: Adds a new command to (un)select random shooting.
///////////////////////////////////////////////////////////////////////////////

#ifndef CCalPrimaryGeneratorMessenger_h
#define CCalPrimaryGeneratorMessenger_h 1

#include "globals.hh"
#include "G4UImessenger.hh"

class CCalPrimaryGeneratorAction;

#include "G4UIcmdWithAString.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"
#include "G4UIcmdWithADouble.hh"
#include "G4UIcmdWithAnInteger.hh"

class CCalPrimaryGeneratorMessenger: public G4UImessenger {
public:
  CCalPrimaryGeneratorMessenger(CCalPrimaryGeneratorAction* myGun);
  ~CCalPrimaryGeneratorMessenger();
  
  void SetNewValue(G4UIcommand * command,G4String newValues);
  
private:
  CCalPrimaryGeneratorAction* myAction;
  
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

