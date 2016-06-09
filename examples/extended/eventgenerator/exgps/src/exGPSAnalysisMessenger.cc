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
/// \file eventgenerator/exgps/src/exGPSAnalysisMessenger.cc
/// \brief Implementation of the exGPSAnalysisMessenger class
//
#ifdef G4ANALYSIS_USE

#include "exGPSAnalysisMessenger.hh"
#include "exGPSAnalysisManager.hh"

#include "G4UIdirectory.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

exGPSAnalysisMessenger::exGPSAnalysisMessenger(exGPSAnalysisManager* analysisManager)
  :fExGPSAnalysis(analysisManager)

{ 
  fExGPSAnalysisDir = new G4UIdirectory("/analysis/");
  fExGPSAnalysisDir->SetGuidance("exGPS analysis control.");
  
  fFileNameCmd = new G4UIcmdWithAString("/analysis/filename",this);
  fFileNameCmd->SetGuidance("Input the name for the AIDA output file.");
  fFileNameCmd->SetParameterName("filename",true,true);
  fFileNameCmd->SetDefaultValue("exgps.aida");

  fFileTypeCmd = new G4UIcmdWithAString("/analysis/filetype",this);
  fFileTypeCmd->SetGuidance("Input the type (xml/hbook/root) for the AIDA output file.");
  fFileTypeCmd->SetParameterName("filetype",true,true);
  fFileTypeCmd->SetDefaultValue("xml");

  fMaxEngCmd = new G4UIcmdWithADoubleAndUnit("/analysis/maxeng",this);
  fMaxEngCmd->SetGuidance("Sets the maximum energy of histo");
  fMaxEngCmd->SetParameterName("maxeng",true,true);
  fMaxEngCmd->SetDefaultUnit("keV");
  fMaxEngCmd->SetUnitCandidates("eV keV MeV GeV TeV PeV");

  fMinEngCmd = new G4UIcmdWithADoubleAndUnit("/analysis/mineng",this);
  fMinEngCmd->SetGuidance("Sets the minimum energy of histo");
  fMinEngCmd->SetParameterName("mineng",true,true);
  fMinEngCmd->SetDefaultUnit("keV");
  fMinEngCmd->SetUnitCandidates("eV keV MeV GeV TeV PeV");

  fMaxPosCmd = new G4UIcmdWithADoubleAndUnit("/analysis/maxpos",this);
  fMaxPosCmd->SetGuidance("Set max length of source position");
  fMaxPosCmd->SetParameterName("maxpos",true,true);
  fMaxPosCmd->SetDefaultUnit("cm");
  fMaxPosCmd->SetUnitCandidates("micron mm cm m km");

  fMinPosCmd = new G4UIcmdWithADoubleAndUnit("/analysis/minpos",this);
  fMinPosCmd->SetGuidance("Set min length of source position");
  fMinPosCmd->SetParameterName("minpos",true,true);
  fMinPosCmd->SetDefaultUnit("cm");
  fMinPosCmd->SetUnitCandidates("micron mm cm m km");

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

exGPSAnalysisMessenger::~exGPSAnalysisMessenger()
{
  delete fFileNameCmd; 
  delete fFileTypeCmd; 
  delete fMaxEngCmd;
  delete fMinEngCmd;
  delete fMaxPosCmd;
  delete fMinPosCmd;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void exGPSAnalysisMessenger::SetNewValue(G4UIcommand* command,
                                                                                                                        G4String newValue)
{ 

  // 1D Histograms

  if( command == fFileNameCmd )
    { 
      fExGPSAnalysis->SetFileName(newValue);
    }
  else if ( command == fFileTypeCmd )
    { 
      fExGPSAnalysis->SetFileType(newValue);
    }
  else if( command == fMaxPosCmd)
    { 
      fExGPSAnalysis->SetPosMax(fMaxPosCmd->GetNewDoubleValue(newValue)); 
    }
  else if( command == fMinPosCmd)
    { 
      fExGPSAnalysis->SetPosMin(fMinPosCmd->GetNewDoubleValue(newValue)); 
    }
  else if( command == fMaxEngCmd)
    { 
      fExGPSAnalysis->SetEngMax(fMaxEngCmd->GetNewDoubleValue(newValue)); 
    }
  else if( command == fMinEngCmd)
    { 
      fExGPSAnalysis->SetEngMin(fMinEngCmd->GetNewDoubleValue(newValue)); 
    }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....


#endif // G4ANALYSIS_USE
