#ifndef __ANALYZER_MESSENGER_DEFINED__
#define __ANALYZER_MESSENGER_DEFINED__

#include "G4UImessenger.hh"

class G4UIdirectory;
class G4UIcmdWithAnInteger;
class G4UIcmdWithADouble;
class G4UIcmdWithoutParameter;
class G4UIparameter;
class G4UIcommand;
class G4UIcmdWithAString;
class Analyzer;

class AnalyzerMessenger : public G4UImessenger
{
public:
  AnalyzerMessenger(Analyzer* pAnalyzer);
  ~AnalyzerMessenger();
  void SetNewValue(G4UIcommand* pCmd,G4String szValue);
private:
  G4UIdirectory*        m_pAnalyzerDir;
  G4UIparameter*        m_pParFileName;
  G4UIparameter*        m_pParIntStyle;
  G4UIparameter*        m_pParHistoID;
  G4UIparameter*        m_pParHistoID1;
  G4UIparameter*        m_pParKey;
  G4UIparameter*        m_pParValue;
  G4UIparameter*        m_pAAOption;
  G4UIparameter*        m_pAAEnergy;
  G4UIparameter*        m_pParStoreID;
  G4UIparameter*        m_pParBias;
  G4UIcmdWithAnInteger* m_pCmdDelete;
  G4UIcommand*          m_pCmdShow;
  G4UIcommand*          m_pCmdRebin;
  G4UIcmdWithAnInteger* m_pCmdHide;
  G4UIcommand*          m_pCmdPS;
  G4UIcmdWithAString*   m_pCmdStore;
  //  G4UIcmdWithADouble*   m_pCmdAA;
  G4UIcommand*          m_pCmdAA;
  G4UIcmdWithADouble*   m_pCmdMinEnergy;
  G4UIcmdWithAnInteger* m_pCmdSave;
  G4UIcmdWithAnInteger* m_pCmdRemove;
  G4UIcommand*          m_pCmdData;
  G4UIcommand*          m_pCmdParameters;
  G4UIcmdWithAnInteger* m_pCmdHistograms;
  G4UIcmdWithAnInteger* m_pCmdBuildDiff;
  G4UIcmdWithoutParameter* m_pCmdReset;
  G4UIcmdWithoutParameter* m_pCmdList;
  G4UIcommand*          m_pCmdReadOptions;
  G4UIcommand*          m_pCmdDeleteOptions;
  G4UIcmdWithAnInteger* m_pCmdListOptions;
  G4UIcmdWithAnInteger* m_pCmdSmooth;
  G4UIcommand*          m_pCmdAdd;

  Analyzer*             m_pAnalyzer;
};
#endif
