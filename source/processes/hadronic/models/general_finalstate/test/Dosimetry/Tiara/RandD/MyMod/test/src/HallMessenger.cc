#include "HallMessenger.hh"
#include "Hall.hh"
#include "G4UIparameter.hh"
#include "G4UIcommand.hh"
#include "G4UIdirectory.hh"
#include "G4UIcmdWithAnInteger.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"
#include "G4UIcmdWithoutParameter.hh"
#include "G4UIcmdWithAString.hh"
#include "G4RunManager.hh"

//#include "stdafx.h"

class RunManager : public G4RunManager
{
public:
  RunManager() : G4RunManager() {;};
  ~RunManager() {;};
  void ProcessTill(G4int nProtons,G4int nNeutrons);
};

extern RunManager* g_pRunManager;

HallMessenger::HallMessenger(Hall* pOwner)
{
  m_pDetector = pOwner;

  m_pDetectorDirectory = new G4UIdirectory("/detector/");
  m_pDetectorDirectory->SetGuidance("resets target and update hall");

  m_pLayerDirectory = new G4UIdirectory("/layer/");
  m_pLayerDirectory->SetGuidance("layer manipulation");

  m_pParLayer = new G4UIparameter("LayerNumber",'i',false);
  m_pParMaterial = new G4UIparameter("LayerMaterial",'s',false);
  m_pParThickness = new G4UIparameter("LayerThickness",'d',false);
  m_pParProtons = new G4UIparameter("ProtonNumber",'i',false);
  m_pParNeutrons = new G4UIparameter("NeutronNumber",'i',false);

  m_pSetLayerCmd =  new G4UIcmdWithAnInteger("/layer/setLayers",this);
  m_pSetLayerCmd->SetParameterName("LayerNumbs",false);
  m_pSetLayerCmd->SetGuidance("Set number of layers in experimental assembly");
  m_pSetLayerCmd->AvailableForStates(Idle);

  m_pCmdDumpImportance = new G4UIcmdWithAString("/layer/print",this);
  m_pCmdDumpImportance->SetGuidance("prints importance and scorer values to a file");
  m_pCmdDumpImportance->SetParameterName("fileName",false);
  m_pCmdDumpImportance->AvailableForStates(Idle);
  m_pImpCmd = new G4UIcmdWithoutParameter("/layer/initImp",this);
  m_pImpCmd->SetGuidance("Initialize importance processes");
  m_pImpCmd->AvailableForStates(Idle);
  m_pNImpCmd = new G4UIcmdWithoutParameter("/layer/delImp",this);
  m_pNImpCmd->SetGuidance("Deletes importance processes");
  m_pNImpCmd->AvailableForStates(Idle);

  m_pSetColimatorCmd = new G4UIcmdWithAnInteger("/layer/setColimator",this);
  m_pSetColimatorCmd->SetGuidance("Sets a layer as a colimator");
  m_pSetColimatorCmd->SetParameterName("ColimatorNum",false);
  m_pSetColimatorCmd->AvailableForStates(Idle);

  m_pSetShieldCmd = new G4UIcmdWithAnInteger("/layer/setShield",this);
  m_pSetShieldCmd->SetGuidance("sets a layer as a shield");
  m_pSetShieldCmd->SetParameterName("ShieldNum",false);
  m_pSetShieldCmd->AvailableForStates(Idle);

  m_pCreateTargetCmd = new G4UIcmdWithADoubleAndUnit("/detector/createTarget",this);
  m_pCreateTargetCmd->SetGuidance("Set target thickness");
  m_pCreateTargetCmd->SetParameterName("TargetThick",true);
  m_pCreateTargetCmd->SetDefaultValue(3.6);
  m_pCreateTargetCmd->SetDefaultUnit("mm");
  m_pCreateTargetCmd->SetUnitCategory("Length");
  m_pCreateTargetCmd->AvailableForStates(Idle);

  m_pUpdateCmd = new G4UIcmdWithoutParameter("/detector/update",this);
  m_pUpdateCmd->SetGuidance("Update detector geometry");
  m_pUpdateCmd->AvailableForStates(Idle);

  m_pResetCmd = new G4UIcmdWithoutParameter("/detector/reset",this);
  m_pResetCmd->SetGuidance("Resets geometry, that is delete layers.");
  m_pResetCmd->AvailableForStates(Idle);

  m_pSetMaterialCmd = new G4UIcommand("/layer/setMat",this);
  m_pSetMaterialCmd->SetParameter(m_pParLayer);
  m_pSetMaterialCmd->SetParameter(m_pParMaterial);
  m_pSetMaterialCmd->SetGuidance("sets material for a layer");
  m_pSetMaterialCmd->AvailableForStates(Idle);

  m_pSetThickCmd = new G4UIcommand("/layer/setThick",this);
  m_pSetThickCmd->SetParameter(m_pParLayer);
  m_pSetThickCmd->SetParameter(m_pParThickness);
  m_pSetThickCmd->AvailableForStates(Idle);

  m_pCmdRunTill = new G4UIcommand("/run/till",this);
  m_pCmdRunTill->SetGuidance("Runs nProtons untill at least nNeutrons presents in all histograms");
  m_pCmdRunTill->SetParameter(m_pParProtons);
  m_pCmdRunTill->SetParameter(m_pParNeutrons);
  m_pCmdRunTill->AvailableForStates(Idle);

  m_pCmdParLayers = new G4UIcmdWithAnInteger("/layer/ParLayers",this);
  m_pCmdParLayers->SetGuidance("Sets n layers in parallel geometry per layer in mass geometry");
  m_pCmdParLayers->SetParameterName("n",true);
  m_pCmdParLayers->SetDefaultValue(1);
  m_pCmdParLayers->AvailableForStates(Idle);
}

HallMessenger::~HallMessenger()
{ 
  //  delete m_pParLayer;
  //  delete m_pParMaterial;
  // delete m_pParThickness;
  delete m_pSetColimatorCmd;
  delete m_pSetShieldCmd;
  delete m_pSetLayerCmd;
  delete m_pCreateTargetCmd;
  delete m_pUpdateCmd;
  delete m_pSetMaterialCmd;
  delete m_pSetThickCmd;
  delete m_pLayerDirectory;
  delete m_pDetectorDirectory;
  delete m_pResetCmd;
  delete m_pCmdRunTill;
  delete m_pCmdDumpImportance;
  delete m_pImpCmd;
  delete m_pNImpCmd;
  delete m_pCmdParLayers;
}
void HallMessenger::SetNewValue(G4UIcommand* pCmd,G4String szVal)
{
  if(pCmd==m_pSetLayerCmd)
    m_pDetector->CreateLayers(m_pSetLayerCmd->GetNewIntValue(szVal));
  else if(pCmd==m_pCreateTargetCmd)
    m_pDetector->CreateTarget(m_pCreateTargetCmd->GetNewDoubleValue(szVal));
  else if(pCmd==m_pUpdateCmd)
    m_pDetector->Update();
  else if(pCmd==m_pResetCmd)
    m_pDetector->ResetGeometry();
  else if(pCmd==m_pSetMaterialCmd){
    G4String szNumber;
    szVal.strip(2);
    G4int pos = szVal.find(' ');
    szNumber = szVal;
    szNumber.remove(pos,szVal.length()-pos);
    szVal.remove(0,pos+1);
    szVal.strip(2);
    pos = atoi(szNumber.data());
    m_pDetector->LayerMat(szVal,pos);
  }
  else if(pCmd == m_pCmdRunTill){
    szVal.strip(2);
    G4String szProtons;
    G4int nProtons,nNeutrons,pos;
    pos = szVal.find(' ');
    szProtons = szVal;
    szProtons.remove(pos,szVal.length()-pos);
    szVal.remove(0,pos);
    nProtons = atoi(szProtons.data());
    nNeutrons = atoi(szVal.data());
    g_pRunManager->ProcessTill(nProtons,nNeutrons);
  }
  else if(pCmd==m_pSetThickCmd){
    szVal.strip(2);
    G4String Number = szVal;
    G4int pos = szVal.find(' ');
    Number.remove(pos,szVal.length()-pos);
    szVal.remove(0,pos);
    szVal.strip(2);
    pos = atoi(Number.data());
    G4double thick = atof(szVal.data());
    m_pDetector->LayerDepth(thick,pos);
  }
  else if(pCmd==m_pSetColimatorCmd)
    m_pDetector->SetAsColimator(m_pSetColimatorCmd->GetNewIntValue(szVal));
  else if(pCmd==m_pSetShieldCmd)
    m_pDetector->SetAsShield(m_pSetShieldCmd->GetNewIntValue(szVal));
  else if(pCmd==m_pCmdDumpImportance){
    m_pDetector->DumpImportance(szVal);
  }
  else if(pCmd==m_pImpCmd){
    m_pDetector->UseImportance();
  }
  else if(pCmd==m_pNImpCmd){
    m_pDetector->DontUseImportance();
  }
  else if(pCmd==m_pCmdParLayers){
    m_pDetector->SetBlocksPerLayer(m_pCmdParLayers->GetNewIntValue(szVal));
  }
}
