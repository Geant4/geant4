#ifndef MLAnalysisMessenger_h
#define MLAnalysisMessenger_h 1
////////////////////////////////////////////////////////////////////////////////
//
#include "globals.hh"
#include "G4UImessenger.hh"

class MLAnalysisManager;
class G4UIdirectory;
class G4UIcmdWithAString;
class G4UIcmdWithAnInteger;
class G4UIcmdWithoutParameter;
class G4UIcmdWithADoubleAndUnit;
class G4UIcmdWithADouble;
////////////////////////////////////////////////////////////////////////////////
//
class MLAnalysisMessenger: public G4UImessenger
{
public:
  MLAnalysisMessenger(MLAnalysisManager* );
  ~MLAnalysisMessenger();
  
  void SetNewValue (G4UIcommand*, G4String);

private:
  MLAnalysisManager          *analysisManager;
  G4UIdirectory              *MLAnalysisDir;
  G4UIcmdWithAString         *FileNameCmd;
  G4UIcmdWithADoubleAndUnit  *NormalCmd;

  G4UIdirectory              *MLAnaFluxDir;
  
  G4UIdirectory              *MLAnaFluEngDir;
  G4UIcmdWithAString         *ModeEngCmd;
  G4UIcmdWithADoubleAndUnit  *MaxEngCmd;
  G4UIcmdWithADoubleAndUnit  *MinEngCmd;
  G4UIcmdWithAnInteger       *NbinEngCmd;
  G4UIcmdWithADoubleAndUnit  *AddEngCmd;
  G4UIcmdWithADoubleAndUnit  *DeleteEngCmd;
  G4UIcmdWithoutParameter    *ListEngCmd;
  G4UIcmdWithoutParameter    *ClearEngCmd;
  G4UIcmdWithoutParameter    *DefaultEngCmd;

  G4UIdirectory              *MLAnaFluAngDir;
  G4UIcmdWithAString         *FluenceTypeCmd;
  G4UIcmdWithAString         *FluenceUnitCmd;
  G4UIcmdWithAString         *ModeAngCmd;
  G4UIcmdWithADoubleAndUnit  *MaxAngCmd;
  G4UIcmdWithADoubleAndUnit  *MinAngCmd;
  G4UIcmdWithAnInteger       *NbinAngCmd;
  G4UIcmdWithADoubleAndUnit  *AddAngCmd;
  G4UIcmdWithADoubleAndUnit  *DeleteAngCmd;
  G4UIcmdWithoutParameter    *ListAngCmd;
  G4UIcmdWithoutParameter    *ClearAngCmd;
  G4UIcmdWithoutParameter    *DefaultAngCmd;

  G4UIdirectory              *MLAnaFluPartDir;
  G4UIcmdWithAString         *AddPartCmd;
  G4UIcmdWithAString         *DeletePartCmd;
  G4UIcmdWithoutParameter    *ListPartCmd;

  G4UIdirectory              *MLAnaPHSDir;
  G4UIdirectory              *MLAnaPHSEngDir;
  G4UIcmdWithAString         *ModePHSCmd;
  G4UIcmdWithADoubleAndUnit  *MaxPHSCmd;
  G4UIcmdWithADoubleAndUnit  *MinPHSCmd;
  G4UIcmdWithAnInteger       *NbinPHSCmd;
  G4UIcmdWithADoubleAndUnit  *AddPHSCmd;
  G4UIcmdWithADoubleAndUnit  *DeletePHSCmd;
  G4UIcmdWithoutParameter    *ListPHSCmd;
  G4UIcmdWithoutParameter    *ClearPHSCmd;
  G4UIcmdWithoutParameter    *DefaultPHSCmd;

  G4UIdirectory              *MLAnaDoseDir;
  G4UIcmdWithAnInteger       *AddDoseCmd;
  G4UIcmdWithAnInteger       *DeleteDoseCmd;
  G4UIcmdWithoutParameter    *ListDoseCmd;
  G4UIcmdWithAString         *DoseUnitCmd;

  G4UIdirectory              *MLAnaNIELDir;
  G4UIcmdWithAnInteger       *AddNIELCmd;
  G4UIcmdWithAnInteger       *DeleteNIELCmd;
  G4UIcmdWithoutParameter    *ListNIELCmd;
  G4UIcmdWithAString         *NIELUnitCmd;
  G4UIcmdWithAString         *SelectNIELCmd;
};
////////////////////////////////////////////////////////////////////////////////
#endif
