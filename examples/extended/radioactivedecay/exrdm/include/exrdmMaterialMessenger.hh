#ifndef exrdmMaterialMessenger_h
#define exrdmMaterialMessenger_h 1
////////////////////////////////////////////////////////////////////////////////
//
#include "globals.hh"
#include "G4UImessenger.hh"

#include "G4UIcommand.hh"
#include "G4UIdirectory.hh"
#include "G4UIparameter.hh"
#include "G4UIcmdWithAnInteger.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIcmdWithoutParameter.hh"

class exrdmMaterial;
////////////////////////////////////////////////////////////////////////////////
//
class exrdmMaterialMessenger: public G4UImessenger
{
public:
  exrdmMaterialMessenger(exrdmMaterial* );
  ~exrdmMaterialMessenger();

  void SetNewValue (G4UIcommand*, G4String);

private:

  exrdmMaterial                *materialsManager;

  G4UIdirectory             *MaterialDir;
  G4UIcmdWithoutParameter   *ListCmd;
  G4UIcmdWithAnInteger      *DeleteIntCmd;
  G4UIcmdWithAString        *DeleteNameCmd;
  G4UIcommand               *AddCmd;
};
////////////////////////////////////////////////////////////////////////////////
#endif
