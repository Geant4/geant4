#ifndef MLMaterialMessenger_h
#define MLMaterialMessenger_h 1
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

class MLMaterial;
////////////////////////////////////////////////////////////////////////////////
//
class MLMaterialMessenger: public G4UImessenger
{
public:
  MLMaterialMessenger(MLMaterial* );
  ~MLMaterialMessenger();

  void SetNewValue (G4UIcommand*, G4String);

private:

  MLMaterial                *materialsManager;

  G4UIdirectory             *MaterialDir;
  G4UIcmdWithoutParameter   *ListCmd;
  G4UIcmdWithAnInteger      *DeleteIntCmd;
  G4UIcmdWithAString        *DeleteNameCmd;
  G4UIcommand               *AddCmd;
};
////////////////////////////////////////////////////////////////////////////////
#endif
