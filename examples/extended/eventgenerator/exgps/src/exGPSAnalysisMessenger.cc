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

#include <AIDA/AIDA.h>

#include "exGPSAnalysisMessenger.hh"
#include "exGPSAnalysisManager.hh"

#include "G4UIdirectory.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

exGPSAnalysisMessenger::exGPSAnalysisMessenger(exGPSAnalysisManager* analysisManager)
  :exGPSAnalysis(analysisManager)

{ 
  exGPSAnalysisDir = new G4UIdirectory("/analysis/");
  exGPSAnalysisDir->SetGuidance("exGPS analysis control.");
  
  FileNameCmd = new G4UIcmdWithAString("/analysis/filename",this);
  FileNameCmd->SetGuidance("Input the name for the AIDA output file.");
  FileNameCmd->SetParameterName("filename",true,true);
  FileNameCmd->SetDefaultValue("exgps.aida");

  FileTypeCmd = new G4UIcmdWithAString("/analysis/filetype",this);
  FileTypeCmd->SetGuidance("Input the type (xml/hbook) for the AIDA output file.");
  FileTypeCmd->SetParameterName("filetype",true,true);
  FileTypeCmd->SetDefaultValue("xml");

  MaxEngCmd = new G4UIcmdWithADoubleAndUnit("/analysis/maxeng",this);
  MaxEngCmd->SetGuidance("Sets the maximum energy of histo");
  MaxEngCmd->SetParameterName("maxeng",true,true);
  MaxEngCmd->SetDefaultUnit("keV");
  MaxEngCmd->SetUnitCandidates("eV keV MeV GeV TeV PeV");

  MinEngCmd = new G4UIcmdWithADoubleAndUnit("/analysis/mineng",this);
  MinEngCmd->SetGuidance("Sets the minimum energy of histo");
  MinEngCmd->SetParameterName("mineng",true,true);
  MinEngCmd->SetDefaultUnit("keV");
  MinEngCmd->SetUnitCandidates("eV keV MeV GeV TeV PeV");

  MaxPosCmd = new G4UIcmdWithADoubleAndUnit("/analysis/maxpos",this);
  MaxPosCmd->SetGuidance("Set max length of source position");
  MaxPosCmd->SetParameterName("maxpos",true,true);
  MaxPosCmd->SetDefaultUnit("cm");
  MaxPosCmd->SetUnitCandidates("micron mm cm m km");

  MinPosCmd = new G4UIcmdWithADoubleAndUnit("/analysis/minpos",this);
  MinPosCmd->SetGuidance("Set min length of source position");
  MinPosCmd->SetParameterName("minpos",true,true);
  MinPosCmd->SetDefaultUnit("cm");
  MinPosCmd->SetUnitCandidates("micron mm cm m km");

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

exGPSAnalysisMessenger::~exGPSAnalysisMessenger()
{
  delete FileNameCmd; 
  delete FileTypeCmd; 
  delete MaxEngCmd;
  delete MinEngCmd;
  delete MaxPosCmd;
  delete MinPosCmd;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void exGPSAnalysisMessenger::SetNewValue(G4UIcommand* command,G4String newValue)
{ 

  // 1D Histograms

  if( command == FileNameCmd )
    { 
      exGPSAnalysis->SetFileName(newValue);
    }
  else if ( command == FileTypeCmd )
    { 
      exGPSAnalysis->SetFileType(newValue);
    }
  else if( command == MaxPosCmd)
    { 
      exGPSAnalysis->SetPosMax(MaxPosCmd->GetNewDoubleValue(newValue)); 
    }
  else if( command == MinPosCmd)
    { 
      exGPSAnalysis->SetPosMin(MinPosCmd->GetNewDoubleValue(newValue)); 
    }
  else if( command == MaxEngCmd)
    { 
      exGPSAnalysis->SetEngMax(MaxEngCmd->GetNewDoubleValue(newValue)); 
    }
  else if( command == MinEngCmd)
    { 
      exGPSAnalysis->SetEngMin(MinEngCmd->GetNewDoubleValue(newValue)); 
    }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....


#endif // G4ANALYSIS_USE






