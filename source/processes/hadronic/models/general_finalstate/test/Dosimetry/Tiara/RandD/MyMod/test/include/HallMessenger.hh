#ifndef __HALL_MESSENGER_DEFINED__
#define __HALL_MESSENGER_DEFINED__

#include "G4UImessenger.hh"
#include "globals.hh"

class G4UIcmdWithAnInteger;
class G4UIcmdWithADoubleAndUnit;
class G4UIparameter;
class G4UIcommand;
class G4UIcmdWithoutParameter;
class Hall;
class G4UIdirectory;
class G4UIcmdWithAString;

class HallMessenger : public G4UImessenger
{
public:
  HallMessenger(Hall* pOwner);
  ~HallMessenger();
  void SetNewValue(G4UIcommand* pCmd,G4String szValue);
private:
  G4UIparameter* m_pParLayer;
  G4UIparameter* m_pParMaterial;
  G4UIparameter* m_pParThickness;
  G4UIparameter* m_pParProtons;
  G4UIparameter* m_pParNeutrons;
  G4UIcmdWithAnInteger*      m_pSetLayerCmd;
  G4UIcmdWithAnInteger*      m_pSetColimatorCmd;
  G4UIcmdWithAnInteger*      m_pSetShieldCmd;
  G4UIcmdWithADoubleAndUnit* m_pCreateTargetCmd;
  G4UIcmdWithoutParameter*   m_pUpdateCmd;
  G4UIcmdWithoutParameter*   m_pResetCmd;
  G4UIcmdWithoutParameter*   m_pImpCmd;
  G4UIcmdWithoutParameter*   m_pNImpCmd;
  G4UIcommand*               m_pSetMaterialCmd;
  G4UIcommand*               m_pSetThickCmd;
  G4UIcommand*               m_pCmdRunTill;
  G4UIcmdWithAString*        m_pCmdDumpImportance;
  G4UIcmdWithAnInteger*      m_pCmdParLayers;
  G4UIdirectory* m_pLayerDirectory;
  G4UIdirectory* m_pDetectorDirectory;
  Hall*          m_pDetector;
};
#endif
