#ifndef MLPhysicsListMessenger_h
#define MLPhysicsListMessenger_h 1
////////////////////////////////////////////////////////////////////////////////
//
#include "globals.hh"
#include "G4UImessenger.hh"

#include "MLPhysicsList.hh"
#include "G4UIdirectory.hh"
#include "G4UIcommand.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIcmdWithoutParameter.hh"
#include "G4UIcmdWithAnInteger.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"
////////////////////////////////////////////////////////////////////////////////
//
class MLPhysicsListMessenger: public G4UImessenger
{
public:
  
  MLPhysicsListMessenger (MLPhysicsList* );
  ~MLPhysicsListMessenger ();
  
  void SetNewValue (G4UIcommand*, G4String);
  
private:
  
  MLPhysicsList* pPhysicsList;
  
  G4UIdirectory*              PhysDir;
  
  G4UIcmdWithAString*         scenarioCmd;
  G4UIcmdWithoutParameter*    listscenCmd;

  G4UIdirectory*              RegionDir;

  G4UIcmdWithAString*         addRegCmd;
  G4UIcmdWithAString*         delRegCmd;
  G4UIcmdWithoutParameter*    listRegCmd;
  G4UIcommand*                addRegLayerCmd;  
  G4UIcommand*                delRegLayerCmd;  
  G4UIcmdWithAString*         listRegLayerCmd;

  G4UIdirectory*              CutsDir;
  G4UIdirectory*              CutsGlobalDir;
  G4UIdirectory*              CutsRegionDir;

  G4UIcmdWithADoubleAndUnit  *gDefaultCmd;
  G4UIcmdWithADoubleAndUnit  *gGammaCmd;
  G4UIcmdWithADoubleAndUnit  *gElectronCmd;
  G4UIcmdWithADoubleAndUnit  *gPositronCmd;

  G4UIcommand*                mSetCutCmd;
  G4UIcommand*                mSetByPartiCmd;

  G4UIcmdWithAnInteger*       verboseCmd;

};
////////////////////////////////////////////////////////////////////////////////
#endif
