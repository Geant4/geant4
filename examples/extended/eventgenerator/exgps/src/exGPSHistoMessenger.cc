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
// $Id: exGPSHistoMessenger.cc 70972 2013-06-07 16:12:12Z gcosmo $
//
/// \file eventgenerator/exgps/src/exGPSHistoMessenger.cc
/// \brief Implementation of the exGPSHistoMessenger class
//



#include "exGPSHistoMessenger.hh"
#include "exGPSHistoManager.hh"

#include "G4UIdirectory.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

exGPSHistoMessenger::exGPSHistoMessenger
                (exGPSHistoManager* histoManager)
  :G4UImessenger(),
   fExGPSHisto(histoManager),
   fFileNameCmd(0),
   fMaxEngCmd(0),
   fMinEngCmd(0),
   fMaxPosCmd(0),
   fMinPosCmd(0)

{ 
  fExGPSHistoDir = new G4UIdirectory("/analysis/");
  fExGPSHistoDir->SetGuidance("exGPS analysis control.");
  
  fFileNameCmd = new G4UIcmdWithAString("/analysis/filename",this);
  fFileNameCmd->SetGuidance("Input the name for the ROOT output file.");
  fFileNameCmd->SetParameterName("filename",true,true);
  fFileNameCmd->SetDefaultValue("exgps");


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

exGPSHistoMessenger::~exGPSHistoMessenger()
{
  delete fFileNameCmd; 
  delete fMaxEngCmd;
  delete fMinEngCmd;
  delete fMaxPosCmd;
  delete fMinPosCmd;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void exGPSHistoMessenger::SetNewValue(G4UIcommand* command,
                                                          G4String newValue)
{ 

  // 1D Histograms

  if( command == fFileNameCmd )
    { 
      fExGPSHisto->SetFileName(newValue);
    }
  else if( command == fMaxPosCmd)
    { 
      fExGPSHisto->SetPosMax(fMaxPosCmd->GetNewDoubleValue(newValue));
    }
  else if( command == fMinPosCmd)
    { 
      fExGPSHisto->SetPosMin(fMinPosCmd->GetNewDoubleValue(newValue));
    }
  else if( command == fMaxEngCmd)
    { 
      fExGPSHisto->SetEMax(fMaxEngCmd->GetNewDoubleValue(newValue));
    }
  else if( command == fMinEngCmd)
    { 
      fExGPSHisto->SetEMin(fMinEngCmd->GetNewDoubleValue(newValue));
    }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......


