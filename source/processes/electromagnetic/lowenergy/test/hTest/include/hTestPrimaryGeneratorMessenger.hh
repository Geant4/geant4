#ifndef hTestPrimaryGeneratorAction_h
#define hTestPrimaryGeneratorAction_h 1

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


// HARP includes
#include "Simulation/HsManager.h"

// G4 includes
#include "globals.hh"
#include "G4UImessenger.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"
#include "G4UIcmdWithAnInteger.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

class HsPrimaryGeneratorOld;

class hTestPrimaryGeneratorAction: public G4UImessenger
{
  public:
  
    hTestPrimaryGeneratorAction(HsPrimaryGeneratorOld* gen);
   ~hTestPrimaryGeneratorAction();
    
    void SetNewValue(G4UIcommand* command, G4String newValue);

  private:
  
    HsPrimaryGeneratorOld*  theGeneratorOld;
    HsManager* theMCmanager;

    G4UIcmdWithADoubleAndUnit* beamXCmd;
    G4UIcmdWithADoubleAndUnit* beamYCmd;
    G4UIcmdWithADoubleAndUnit* beamZCmd;
    G4UIcmdWithADoubleAndUnit* sigmaXCmd;
    G4UIcmdWithADoubleAndUnit* sigmaYCmd;
    G4UIcmdWithADoubleAndUnit* sigmaZCmd;
    G4UIcmdWithADoubleAndUnit* sigmaECmd;
    G4UIcmdWithADoubleAndUnit* maxThetaCmd;

};

#endif

