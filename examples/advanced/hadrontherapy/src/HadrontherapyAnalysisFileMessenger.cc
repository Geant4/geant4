

#include "HadrontherapyAnalysisFileMessenger.hh"
#include "HadrontherapyAnalysisManager.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIcmdWithABool.hh"
#include "G4UIdirectory.hh"

#include "HadrontherapyMatrix.hh"

//
/////////////////////////////////////////////////////////////////////////////
HadrontherapyAnalysisFileMessenger::HadrontherapyAnalysisFileMessenger(HadrontherapyAnalysisManager* amgr)
:AnalysisManager(amgr)
{  

  secondariesCmd = new G4UIcmdWithABool("/analysis/secondaries",this);
  secondariesCmd -> SetParameterName("secondaries", false);
  secondariesCmd -> SetGuidance("If true dose/fluence for the secondaries are written"); 
  secondariesCmd -> AvailableForStates(G4State_Idle,G4State_PreInit);
    // With this messenger you can:
    // give a name to the generated .root file
    // close the old .root file and create a new one (with the given name)
    // One can use this messenger to define a different .root file name other thaen the default one 
#ifdef G4ANALYSIS_USE_ROOT
  FileNameCmd = new G4UIcmdWithAString("/analysis/setAnalysisFile",this);
  FileNameCmd->SetGuidance("Set the .root filename for the root-output");
  FileNameCmd->SetDefaultValue("default.root");
  FileNameCmd->SetParameterName("choice",true); ///<doc did not say what second boolean really does
  FileNameCmd->AvailableForStates(G4State_Idle,G4State_PreInit);
#endif
}

/////////////////////////////////////////////////////////////////////////////
HadrontherapyAnalysisFileMessenger::~HadrontherapyAnalysisFileMessenger()
{
  delete secondariesCmd; 
#ifdef G4ANALYSIS_USE_ROOT
  delete FileNameCmd;
#endif
}

/////////////////////////////////////////////////////////////////////////////
void HadrontherapyAnalysisFileMessenger::SetNewValue(G4UIcommand* command, G4String newValue)
{ 
    if (command == secondariesCmd)
    {
      if (HadrontherapyMatrix::GetInstance())
	  HadrontherapyMatrix::GetInstance() -> secondaries = secondariesCmd -> GetNewBoolValue(newValue);
    }
#ifdef G4ANALYSIS_USE_ROOT
    if (command == FileNameCmd)
    {
	AnalysisManager->SetAnalysisFileName(newValue);
    }
#endif
}

