//
// ********************************************************************
// * DISCLAIMER                                                       *
// *                                                                  *
// * The following disclaimer summarizes all the specific disclaimers *
// * of contributors to this software. The specific disclaimers,which *
// * govern, are listed with their locations in:                      *
// *   http://cern.ch/geant4/license                                  *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.                                                             *
// *                                                                  *
// * This  code  implementation is the  intellectual property  of the *
// * GEANT4 collaboration.                                            *
// * By copying,  distributing  or modifying the Program (or any work *
// * based  on  the Program)  you indicate  your  acceptance of  this *
// * statement, and all its terms.                                    *
// ********************************************************************
//
#ifdef G4ANALYSIS_USE
#include "fluoTestAnalysisMessenger.hh"

#include "fluoTestAnalysisManager.hh"
#include "G4UIdirectory.hh"
#include "G4UIcmdWithAString.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

fluoTestAnalysisMessenger::fluoTestAnalysisMessenger(fluoTestAnalysisManager* analysisManager)
  :fluoTestAnalysis(analysisManager)

{ 
  fluoTestAnalysisDir = new G4UIdirectory("/analysis/");
  fluoTestAnalysisDir->SetGuidance("esperimento analysis control.");
  
     
  //   The Draw command gives the possibility to draw the 1d histograms 
  //   at every event.

  //   The Save command gives the possibility to save the 1d histograms in
  //   two separate PostScript files at the end of the run.
 

  Histo1DDrawCmd = new G4UIcmdWithAString("/analysis/histo1dDraw",this);
  Histo1DDrawCmd->SetGuidance("Enable the drawing of the 1d histograms every event.");
  Histo1DDrawCmd->SetGuidance("Choice: disable, enable(default)");
  Histo1DDrawCmd->SetParameterName("choice",true);
  Histo1DDrawCmd->SetDefaultValue("enable");
  Histo1DDrawCmd->SetCandidates("disable enable");
  Histo1DDrawCmd->AvailableForStates(Idle);

  Histo1DSaveCmd = new G4UIcmdWithAString("/analysis/histo1dSave",this);
  Histo1DSaveCmd->SetGuidance("Enable the saving of the 1d histograms every run.");
  Histo1DSaveCmd->SetGuidance("Choice: disable, enable(default)");
  Histo1DSaveCmd->SetParameterName("choice",true);
  Histo1DSaveCmd->SetDefaultValue("enable");
  Histo1DSaveCmd->SetCandidates("disable enable");
  Histo1DSaveCmd->AvailableForStates(Idle);

 
  }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

fluoTestAnalysisMessenger::~fluoTestAnalysisMessenger()
{
  delete Histo1DDrawCmd; 
  delete Histo1DSaveCmd; 
 
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void fluoTestAnalysisMessenger::SetNewValue(G4UIcommand* command,G4String newValue)
{ 

  // 1D Histograms

  if( command == Histo1DDrawCmd )
    { fluoTestAnalysis->SetHisto1DDraw(newValue);}

  if( command == Histo1DSaveCmd )
    { fluoTestAnalysis->SetHisto1DSave(newValue);}

 
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#endif













