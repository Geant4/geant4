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
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "exrdmHistoMessenger.hh"

#include <sstream>

#include "exrdmHisto.hh"
#include "G4UIdirectory.hh"
#include "G4UIcommand.hh"
#include "G4UIparameter.hh"
#include "G4UIcmdWithAString.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

exrdmHistoMessenger::exrdmHistoMessenger(exrdmHisto* manager)
:histo (manager)
{
  histoDir = new G4UIdirectory("/histo/");
  histoDir->SetGuidance("histograms control");

  factoryCmd = new G4UIcmdWithAString("/histo/fileName",this);
  factoryCmd->SetGuidance("set name for the histograms file");
  factoryCmd->AvailableForStates(G4State_PreInit,G4State_Idle); 

  fileCmd = new G4UIcmdWithAString("/histo/fileType",this);
  fileCmd->SetGuidance("set type (aida, hbook or root) for the histograms file");
  fileCmd->AvailableForStates(G4State_PreInit,G4State_Idle); 

  histoCmd = new G4UIcommand("/histo/setHisto",this);
  histoCmd->SetGuidance("Set bining of the histo number ih :");
  histoCmd->SetGuidance("  nbBins; valMin; valMax; unit (of vmin and vmax)");
  histoCmd->AvailableForStates(G4State_PreInit,G4State_Idle); 

  //
  G4UIparameter* ih = new G4UIparameter("ih",'i',false);
  ih->SetGuidance("histo number : from 0 to MaxexrdmHisto");
  ih->SetParameterRange("ih>=0");
  histoCmd->SetParameter(ih);
  //
  G4UIparameter* nbBins = new G4UIparameter("nbBins",'i',false);
  nbBins->SetGuidance("number of bins");
  nbBins->SetParameterRange("nbBins>0");
  histoCmd->SetParameter(nbBins);
  //
  G4UIparameter* valMin = new G4UIparameter("valMin",'d',false);
  valMin->SetGuidance("valMin, expressed in unit");
  histoCmd->SetParameter(valMin);
  //
  G4UIparameter* valMax = new G4UIparameter("valMax",'d',false);
  valMax->SetGuidance("valMax, expressed in unit");
  histoCmd->SetParameter(valMax);
  //
  G4UIparameter* unit = new G4UIparameter("unit",'s',true);
  unit->SetGuidance("if omitted, vmin and vmax are assumed dimensionless");
  unit->SetDefaultValue("none");
  histoCmd->SetParameter(unit);

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

exrdmHistoMessenger::~exrdmHistoMessenger()
{
  delete fileCmd;
  delete histoCmd;
  delete factoryCmd;
  delete histoDir;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void exrdmHistoMessenger::SetNewValue(G4UIcommand* command,G4String newValues)
{
  if (command == factoryCmd)
    histo->setFileName(newValues);

  if (command == fileCmd)
    histo->setFileType(newValues);
    
  if (command == histoCmd)
   { G4int ih,nbBins; G4double vmin,vmax; char unts[30];
     const char* t = newValues;
     std::istringstream is(t);
     is >> ih >> nbBins >> vmin >> vmax >> unts;
     G4String unit = unts;
     G4double vUnit = 1. ;
     if (unit != "none") vUnit = G4UIcommand::ValueOf(unit);
     histo->setHisto1D(ih,nbBins,vmin,vmax,vUnit);
   }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
