//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
// HadrontherapyAnalysisFileMessenger.cc
//
// See more at: http://g4advancedexamples.lngs.infn.it/Examples/hadrontherapy
//


#include "HadrontherapyAnalysisFileMessenger.hh"
#include "HadrontherapyAnalysisManager.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIcmdWithABool.hh"
#include "G4UIdirectory.hh"

#include "HadrontherapyMatrix.hh"

/////////////////////////////////////////////////////////////////////////////
HadrontherapyAnalysisFileMessenger::HadrontherapyAnalysisFileMessenger(HadrontherapyAnalysisManager* amgr)
:AnalysisManager(amgr)
{  
  secondariesCmd = new G4UIcmdWithABool("/analysis/secondaries",this);
  secondariesCmd -> SetParameterName("secondaries", true);
  secondariesCmd -> SetDefaultValue("true");
  secondariesCmd -> SetGuidance("Set if dose/fluence for the secondaries are written"); 
  secondariesCmd -> AvailableForStates(G4State_Idle, G4State_PreInit);

  // With this messenger you can:
  // give a name to the generated .root file
  // One can use this messenger to define a different .root file name other then the default one 
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

