#ifndef hTestPrimaryGeneratorMessenger_h
#define hTestPrimaryGeneratorMessenger_h 1

//---------------------------------------------------------------------------
//
// ClassName:   hTestPrimaryGeneratorAction
//  
// Description: Definition of physics list parameters via UI interface
//
// Author:      V.Ivanchenko 26/09/00
//
//----------------------------------------------------------------------------
//

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#include "globals.hh"
#include "G4UImessenger.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"
#include "G4UIcmdWithAnInteger.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

class hTestPrimaryGeneratorAction;

class hTestPrimaryGeneratorMessenger: public G4UImessenger
{
  public:
  
    hTestPrimaryGeneratorMessenger(hTestPrimaryGeneratorAction* gen);
   ~hTestPrimaryGeneratorMessenger();
    
    void SetNewValue(G4UIcommand* command, G4String newValue);

  private:
  
    hTestPrimaryGeneratorAction*  theGen;

    G4UIcmdWithADoubleAndUnit* beamXCmd;
    G4UIcmdWithADoubleAndUnit* beamYCmd;
    G4UIcmdWithADoubleAndUnit* beamZCmd;
    G4UIcmdWithADoubleAndUnit* beamECmd;
    G4UIcmdWithADoubleAndUnit* sigmaXCmd;
    G4UIcmdWithADoubleAndUnit* sigmaYCmd;
    G4UIcmdWithADoubleAndUnit* sigmaZCmd;
    G4UIcmdWithADoubleAndUnit* sigmaECmd;
    G4UIcmdWithADoubleAndUnit* maxThetaCmd;

};

#endif

