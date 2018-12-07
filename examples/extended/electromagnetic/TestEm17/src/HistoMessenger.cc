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
/// \file electromagnetic/TestEm17/src/HistoMessenger.cc
/// \brief Implementation of the HistoMessenger class
//
//
// 
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo...... 

#include "HistoMessenger.hh"

#include <sstream>

#include "HistoManager.hh"
#include "G4UIdirectory.hh"
#include "G4UIcommand.hh"
#include "G4UIparameter.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIcmdWithAnInteger.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

HistoMessenger::HistoMessenger(HistoManager* manager)
:G4UImessenger(),fHistoManager (manager),
 fHistoDir(0),   
 fFileNameCmd(0),
 fHistoCmd(0),
 fPrtHistoCmd(0)
{
  fHistoDir = new G4UIdirectory("/testem/histo/");
  fHistoDir->SetGuidance("histograms control");

  fFileNameCmd = new G4UIcmdWithAString("/testem/histo/setFileName",this);
  fFileNameCmd->SetGuidance("set name for the histograms file");
  
  fHistoCmd = new G4UIcommand("/testem/histo/setHisto",this);
  fHistoCmd->SetGuidance("Set bining of the histo number ih :");
  fHistoCmd->SetGuidance("  nbBins; valMin; valMax; unit (of vmin and vmax)");
  //
  G4UIparameter* ih = new G4UIparameter("ih",'i',false);
  ih->SetGuidance("histo number : from 1 to kMaxHisto");
  ih->SetParameterRange("ih>0");
  fHistoCmd->SetParameter(ih);
  //
  G4UIparameter* nbBins = new G4UIparameter("nbBins",'i',false);
  nbBins->SetGuidance("number of bins");
  nbBins->SetParameterRange("nbBins>0");
  fHistoCmd->SetParameter(nbBins);
  //
  G4UIparameter* valMin = new G4UIparameter("valMin",'d',false);
  valMin->SetGuidance("valMin, expressed in unit");
  fHistoCmd->SetParameter(valMin);  
  //    
  G4UIparameter* valMax = new G4UIparameter("valMax",'d',false);
  valMax->SetGuidance("valMax, expressed in unit");
  fHistoCmd->SetParameter(valMax);
  //    
  G4UIparameter* unit = new G4UIparameter("unit",'s',true);
  unit->SetGuidance("if omitted, vmin and vmax are assumed dimensionless");
  unit->SetDefaultValue("none");
  fHistoCmd->SetParameter(unit);
  
  fPrtHistoCmd = new G4UIcmdWithAnInteger("/testem/histo/printHisto",this);
  fPrtHistoCmd->SetGuidance("print histo #id on ascii file");
  fPrtHistoCmd->SetParameterName("id",false);
  fPrtHistoCmd->SetRange("id>0");
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

HistoMessenger::~HistoMessenger()
{
  delete fPrtHistoCmd;  
  delete fHistoCmd;
  delete fFileNameCmd;
  delete fHistoDir;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void HistoMessenger::SetNewValue(G4UIcommand* command, G4String newValues)
{
  if (command == fFileNameCmd)
    fHistoManager->SetFileName(newValues);

  if (command == fHistoCmd)
   { G4int ih,nbBins; G4double vmin,vmax; G4String unit;
     const char* t = newValues;
     std::istringstream is(t);
     is >> ih >> nbBins >> vmin >> vmax >> unit;
     fHistoManager->SetHisto (ih,nbBins,vmin,vmax,unit);     
   }
   
  if (command == fPrtHistoCmd)
    fHistoManager->PrintHisto(fPrtHistoCmd->GetNewIntValue(newValues));   
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

