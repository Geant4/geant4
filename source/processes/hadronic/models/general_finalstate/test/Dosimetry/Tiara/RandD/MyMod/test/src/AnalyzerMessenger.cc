#include "AnalyzerMessenger.hh"
#include "G4UIdirectory.hh"
#include "G4UIparameter.hh"
#include "G4UIcommand.hh"
#include "G4UIcmdWithAnInteger.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIcmdWithADouble.hh"
#include "G4UIcmdWithoutParameter.hh"
#include "NeutAnalyzer.hh"
#include "AnalyzerLin.hh"
#include "string.h"

//#include "stdafx.h"
AnalyzerMessenger::AnalyzerMessenger(Analyzer* pAnalyzer)
{
  m_pAnalyzer = pAnalyzer;

  m_pAnalyzerDir = new G4UIdirectory("/histo/");
  m_pAnalyzerDir->SetGuidance("Commands for managing histograms");

  m_pParFileName = new G4UIparameter("FileName",'s',false);
  m_pParStoreID = new G4UIparameter("Stored histo id",'i',false);
  m_pParIntStyle = new G4UIparameter("PlotStyle",'i',true);
  m_pParIntStyle->SetDefaultValue("0");
  m_pParHistoID = new G4UIparameter("INTHISTOID",'i',false);
  m_pParHistoID1 = new G4UIparameter("INTHISTOID1",'i',false);
  m_pParBias = new G4UIparameter("INTBIAS",'i',false);
  m_pParKey = new G4UIparameter("KeyValue",'s',false);
  m_pParValue = new G4UIparameter("DataValue",'s',false);
  m_pAAOption = new G4UIparameter("AAOption",'s',true);
  m_pAAOption->SetDefaultValue("none");
  m_pAAEnergy = new G4UIparameter("AAEnergy",'d',true);
  m_pAAEnergy->SetDefaultValue("100");

  m_pCmdDelete = new G4UIcmdWithAnInteger("/histo/del",this);
  m_pCmdDelete->SetGuidance("delete a histogram by histoID");
  m_pCmdDelete->SetParameterName("HISTOID",false);
  m_pCmdDelete->AvailableForStates(Idle,GeomClosed);

  m_pCmdShow = new G4UIcommand("/histo/show",this);
  m_pCmdShow->SetGuidance("show a histogram by histoID");
  m_pCmdShow->SetParameter(m_pParHistoID);
  m_pCmdShow->SetParameter(m_pParIntStyle);
  m_pCmdShow->AvailableForStates(Idle,GeomClosed);

  m_pCmdRebin = new G4UIcommand("/histo/rebin",this);
  m_pCmdRebin->SetGuidance("rebining histogram according to bins of first additional histogram");
  m_pCmdRebin->SetParameter(m_pParHistoID);
  m_pCmdRebin->SetParameter(m_pParHistoID1);
  m_pCmdRebin->SetParameter(m_pParBias);
  m_pCmdRebin->AvailableForStates(Idle);

  m_pCmdHide = new G4UIcmdWithAnInteger("/histo/hide",this);
  m_pCmdHide->SetParameterName("HISTOID",false);
  m_pCmdHide->SetGuidance("Hide previously shown histogram");
  m_pCmdHide->AvailableForStates(Idle,GeomClosed);

  m_pCmdPS = new G4UIcommand("/histo/ps",this);
  m_pCmdPS->SetGuidance("Sends histogram view to a ps file");
  m_pCmdPS->SetParameter(m_pParHistoID);
  m_pCmdPS->SetParameter(m_pParFileName);

  m_pCmdStore = new G4UIcmdWithAString("/histo/storeName",this);
  m_pCmdStore->SetGuidance("Choose a store name");
  m_pCmdStore->SetParameterName("FILENAME",false);
  m_pCmdStore->AvailableForStates(Idle,GeomClosed);

  //  m_pCmdAA = new G4UIcmdWithADouble("/histo/setAA",this);
  m_pCmdAA = new G4UIcommand("/histo/setAA",this);
  m_pCmdAA->SetGuidance("Trying to smooth histograms before calculating differences.");
  m_pCmdAA->SetGuidance("Energy !=0 - set antialiasing over 2 neightboring bins below that energy (in MeV)");
  m_pCmdAA->SetGuidance("==0 - use not antialiased data");
  m_pCmdAA->SetGuidance("AA types:");
  m_pCmdAA->SetGuidance("   none - does not smoothing");
  m_pCmdAA->SetGuidance("   box - using box filter");
  m_pCmdAA->SetGuidance("   tent - using tent(triangle) filter");
  m_pCmdAA->SetParameter(m_pAAOption);
  m_pCmdAA->SetParameter(m_pAAEnergy);
  m_pCmdAA->AvailableForStates(Idle);

  m_pCmdMinEnergy = new G4UIcmdWithADouble("/histo/setMinEnergy",this);
  m_pCmdMinEnergy->SetGuidance("Set minimum energy when building the differences");
  m_pCmdMinEnergy->SetParameterName("MinEnergy",true);
  m_pCmdMinEnergy->SetDefaultValue(0.);
  m_pCmdMinEnergy->AvailableForStates(Idle);

  m_pCmdSave = new G4UIcmdWithAnInteger("/histo/save",this);
  m_pCmdSave->SetGuidance("Save a histo to the selected store");
  m_pCmdSave->SetParameterName("HISTOID",false);
  m_pCmdSave->AvailableForStates(Idle);

  m_pCmdRemove = new G4UIcmdWithAnInteger("/histo/remove",this);
  m_pCmdRemove->SetGuidance("Remove a histo from selected store");
  m_pCmdRemove->SetParameterName("HISTOID",false);
  m_pCmdRemove->AvailableForStates(Idle,GeomClosed);

  m_pCmdData = new G4UIcommand("/histo/data",this);
  m_pCmdData->SetGuidance("Add dataset to a histogram");
  m_pCmdData->SetParameter(m_pParHistoID);
  m_pCmdData->SetParameter(m_pParFileName);
  m_pCmdData->AvailableForStates(Idle,GeomClosed);

  m_pCmdParameters = new G4UIcommand("/histo/parameter",this);
  m_pCmdParameters->SetGuidance("Set a parameter of a histogram");
  m_pCmdParameters->SetParameter(m_pParHistoID);
  m_pCmdParameters->SetParameter(m_pParKey);
  m_pCmdParameters->SetParameter(m_pParValue);

  m_pCmdHistograms = new G4UIcmdWithAnInteger("/histo/use",this);
  m_pCmdHistograms->SetGuidance("Create histograms");
  m_pCmdHistograms->SetGuidance("     1 - use at beam detector");
  m_pCmdHistograms->SetGuidance("     2 - use detector 20 cm from beam");
  m_pCmdHistograms->SetGuidance("     4 - use detector 40 cm from beam");
  m_pCmdHistograms->SetGuidance("     8 - use ghost detector( distribution from target)");
  m_pCmdHistograms->SetGuidance(" the parameter is bitwise OR from these values :-)");
  m_pCmdHistograms->SetParameterName("HistoEnum",true);
  m_pCmdHistograms->SetRange("HistoEnum > 0 && HistoEnum < 16");
  m_pCmdHistograms->SetDefaultValue(15);
  m_pCmdHistograms->AvailableForStates(Idle);

  m_pCmdBuildDiff = new G4UIcmdWithAnInteger("/histo/diff",this);
  m_pCmdBuildDiff->SetGuidance("Build differences for given histogram if data are present");
  m_pCmdBuildDiff->SetParameterName("HISTOENUM",false);
  m_pCmdBuildDiff->SetRange("HISTOENUM > 0 && HISTOENUM < 16");
  m_pCmdBuildDiff->AvailableForStates(Idle);

  m_pCmdReset = new G4UIcmdWithoutParameter("/histo/reset",this);
  m_pCmdReset->SetGuidance("Deletes all of current histograms");
  m_pCmdReset->AvailableForStates(Idle);

  m_pCmdList = new G4UIcmdWithoutParameter("/histo/list",this);
  m_pCmdList->SetGuidance("Shows histograms stored in memory");
  m_pCmdList->AvailableForStates(Idle);

  m_pCmdReadOptions = new G4UIcommand("/histo/readOptions",this);
  m_pCmdReadOptions->SetGuidance("Reads histogram options from file");
  m_pCmdReadOptions->SetParameter(m_pParHistoID);
  m_pCmdReadOptions->SetParameter(m_pParFileName);
  m_pCmdReadOptions->AvailableForStates(Idle);

  m_pCmdDeleteOptions = new G4UIcommand("/histo/delOption",this);
  m_pCmdDeleteOptions->SetGuidance("Deletes option assigned to a histogram");
  m_pCmdDeleteOptions->SetParameter(m_pParHistoID);
  m_pCmdDeleteOptions->SetParameter(m_pParKey);
  m_pCmdDeleteOptions->AvailableForStates(Idle);

  m_pCmdListOptions = new G4UIcmdWithAnInteger("/histo/options",this);
  m_pCmdListOptions->SetGuidance("Shows options attached to a histo");
  m_pCmdListOptions->SetParameterName("HistoID",false);
  m_pCmdListOptions->AvailableForStates(Idle);

  m_pCmdSmooth = new G4UIcmdWithAnInteger("/histo/smooth",this);
  m_pCmdSmooth->SetGuidance("Smoothing a histogram");
  m_pCmdSmooth->SetParameterName("HistoID",false);
  m_pCmdSmooth->AvailableForStates(Idle);

  m_pCmdAdd = new G4UIcommand("/histo/add",this);
  m_pCmdAdd->SetGuidance("Add histograms to be plot with a created histogram");
  m_pCmdAdd->SetParameter(m_pParHistoID);
  m_pCmdAdd->SetParameter(m_pParFileName);
  m_pCmdAdd->SetParameter(m_pParStoreID);
  m_pCmdAdd->AvailableForStates(Idle);
}

AnalyzerMessenger::~AnalyzerMessenger()
{
  delete m_pCmdDeleteOptions;
  delete m_pCmdListOptions;
  delete m_pCmdDelete;
  delete m_pCmdShow;
  delete m_pCmdHide;
  delete m_pCmdPS;
  delete m_pCmdStore;
  delete m_pCmdAA;
  delete m_pCmdMinEnergy;
  delete m_pCmdSave;
  delete m_pCmdRemove;
  // delete m_pCmdData;
  //delete m_pCmdReadOptions;
  delete m_pCmdParameters;
  delete m_pCmdHistograms;
  delete m_pCmdBuildDiff;
  delete m_pCmdReset;
  delete m_pCmdList;
  //  delete m_pCmdSmooth;
  delete m_pAnalyzerDir;
}

void AnalyzerMessenger::SetNewValue(G4UIcommand* pCmd,G4String szValue)
{
  if(pCmd==m_pCmdReadOptions){
    G4String szHistoID;
    szValue.strip(2);
    G4int pos = szValue.find(' ');
    szHistoID = szValue;
    szHistoID.remove(pos,szValue.length()-pos);
    szValue.remove(0,pos+1);
    pos = atoi(szHistoID.data());
    szValue.strip(2);
    m_pAnalyzer->ReadOptions(pos,szValue);
  }
  else if(pCmd==m_pCmdDeleteOptions){
    G4String szHistoID;
    szValue.strip(2);
    G4int pos = szValue.find(' ');
    szHistoID = szValue;
    szHistoID.remove(pos,szValue.length()-pos);
    szValue.remove(0,pos+1);
    pos = atoi(szHistoID.data());
    szValue.strip(2);
    m_pAnalyzer->DeleteOption(pos,szValue);
  }
  else if(pCmd==m_pCmdListOptions)
    m_pAnalyzer->ListOptions(m_pCmdListOptions->GetNewIntValue(szValue));
  else if(pCmd==m_pCmdList)
    alListHistos();
  else if(pCmd==m_pCmdDelete)
    m_pAnalyzer->DeleteHisto(m_pCmdDelete->GetNewIntValue(szValue));
  else if(pCmd == m_pCmdReset)
    m_pAnalyzer->Reset();
  else if(pCmd==m_pCmdShow){
    G4int histo,style;
    G4String szHistoID;
    szValue.strip(2);
    G4int pos = szValue.find(' ');
    szHistoID = szValue;
    szHistoID.remove(pos,szValue.length()-pos);
    szValue.remove(0,pos);
    szValue.strip(2);
    histo = atoi(szHistoID.data());
    style = atoi(szValue.data());
    m_pAnalyzer->ShowHisto(histo);
  }
  else if(pCmd==m_pCmdHide)
    m_pAnalyzer->HideHisto(m_pCmdHide->GetNewIntValue(szValue));
  else if(pCmd == m_pCmdPS){
    G4String szHistoID;
    szValue.strip(2);
    G4int pos = szValue.find(' ');
    szHistoID = szValue;
    szHistoID.remove(pos,szValue.length()-pos);
    szValue.remove(0,pos+1);
    pos = atoi(szHistoID.data());
    szValue.strip(2);
    m_pAnalyzer->MakePS(pos,szValue);
  }
  else if(pCmd==m_pCmdStore)
    alOpenStorage(szValue);
  else if(pCmd == m_pCmdSave)
    m_pAnalyzer->SaveHisto(m_pCmdSave->GetNewIntValue(szValue));
  else if(pCmd == m_pCmdRemove)
    alRemoveFromStore(m_pCmdRemove->GetNewIntValue(szValue));
  else if(pCmd == m_pCmdData){
    G4String szHistoID;
    szValue.strip(2);
    G4int pos = szValue.find(' ');
    szHistoID = szValue;
    szHistoID.remove(pos,szValue.length()-pos);
    szValue.remove(0,pos+1);
    szValue.strip(2);
    pos = atoi(szHistoID.data());
    m_pAnalyzer->AddData(pos,szValue);
  }
  else if(pCmd==m_pCmdParameters){
    G4String szHistoID;
    G4String szKey;
    szValue.strip(2);
    G4int pos = szValue.find(' ');
    szHistoID = szValue;
    szHistoID.remove(pos,szValue.length()-pos);
    szValue.remove(0,pos+1);
    szValue.strip(2);
    pos = szValue.find(' ');
    szKey = szValue;
    szKey.remove(pos,szValue.length()-pos);
    szValue.remove(0,pos+1);
    szValue.strip(2);
    pos = atoi(szHistoID);
    m_pAnalyzer->SetProperty(pos,szKey,szValue);
  }
  else if(pCmd==m_pCmdHistograms)
    m_pAnalyzer->Histograms(m_pCmdHistograms->GetNewIntValue(szValue));
  else if(pCmd == m_pCmdBuildDiff){
    int value = m_pCmdBuildDiff->GetNewIntValue(szValue);
    if(value & 8)
      m_pAnalyzer->BuildDifferences(8);
    if(value & 4)
      m_pAnalyzer->BuildDifferences(4);
    if(value & 2)
      m_pAnalyzer->BuildDifferences(2);
    if(value & 1)
      m_pAnalyzer->BuildDifferences(1);
  }
  else if(pCmd==m_pCmdAA){
    G4int pos = (szValue.strip(2),szValue.find(' '));
    G4String str = szValue;
    str.remove(pos,szValue.length()-pos);
    szValue.remove(0,pos+1);
    G4double en = atof(szValue.data());
    str.toLower();
    if(str==G4String("none"))
      alSetAA(AA_NONE);
    else if(str==G4String("box"))
      alSetAA(AA_BOX);
    else if(str==G4String("tent"))
      alSetAA(AA_TENT);
    else{
      G4cout<<"error: "<<str<<" is not known"<<G4endl;
      return;
    }
    m_pAnalyzer->UseCut(en);
  }
  else if(pCmd==m_pCmdMinEnergy){
    double dMinEnergy = m_pCmdMinEnergy->GetNewDoubleValue(szValue);
    m_pAnalyzer->SetMinEnergy(dMinEnergy);
  }
  else if(pCmd==m_pCmdSmooth){
    m_pAnalyzer->BuildSmooth(m_pCmdSmooth->GetNewIntValue(szValue));
  }
  else if(pCmd==m_pCmdAdd){
    int HistoID;
    int StoreID;
    int pos;
    szValue.strip(2);
    G4String szTag;
    pos = szValue.find(' ');
    szTag = szValue;
    szTag.remove(pos,szValue.length()-pos);
    szValue.remove(0,pos+1);
    szValue.strip(2);
    HistoID = atoi(szTag.data());
    pos = szValue.find(' ');
    szTag = szValue;
    szTag.remove(pos,szValue.length()-pos);
    szTag.strip(2);
    szValue.remove(0,pos+1);
    szValue.strip(2);
    if(szValue==""){
      StoreID = atoi(szTag.data());
      szTag = szValue;
    }
    else
      StoreID = atoi(szValue.data());
    m_pAnalyzer->AdditionalHisto(HistoID,szTag,StoreID);
  }
  else if(pCmd == m_pCmdRebin){
    G4String szStart = szValue;
    G4String szEnd;
    unsigned start,end,bias,pos;
    szStart.strip(2);     szValue.strip(2);
    pos = szStart.find(' ');
    szStart.remove(pos,szStart.length()-pos);
    szValue.remove(0,pos+1);
    szStart.strip(2);
    szValue.strip(2);
    szEnd = szValue;
    pos = szEnd.find(' ');
    szEnd.remove(pos,szEnd.length()-pos);
    szValue.remove(0,pos+1);
    szEnd.strip(2);
    szValue.strip(2);
    start = atoi(szStart.c_str());
    end = atoi(szEnd.c_str());
    bias = atoi(szValue.c_str());
    m_pAnalyzer->Rebin(start,end,bias);
    
  }
}
