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
//
// $Id: HistoMessenger.cc,v 1.2 2010/09/17 11:21:46 maire Exp $
// GEANT4 tag $Name: geant4-09-04 $
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
:histoManager (manager)
{
  histoDir = new G4UIdirectory("/rdecay01/histo/");
  histoDir->SetGuidance("histograms control");

  factoryCmd = new G4UIcmdWithAString("/rdecay01/histo/setFileName",this);
  factoryCmd->SetGuidance("set name for the histograms file");

  typeCmd = new G4UIcmdWithAString("/rdecay01/histo/setFileType",this);
  typeCmd->SetGuidance("set histograms file type: hbook, root, XML");
  typeCmd->SetCandidates("hbook root XML");

  optionCmd = new G4UIcmdWithAString("/rdecay01/histo/setFileOption",this);
  optionCmd->SetGuidance("set option for the histograms file");
  
  histoCmd = new G4UIcommand("/rdecay01/histo/setHisto",this);
  histoCmd->SetGuidance("Set bining of the histo number ih :");
  histoCmd->SetGuidance("  nbBins; valMin; valMax; unit (of vmin and vmax)");
  //
  G4UIparameter* ih = new G4UIparameter("ih",'i',false);
  ih->SetGuidance("histo number : from 1 to MaxHisto");
  ih->SetParameterRange("ih>0");
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
  
  prhistoCmd = new G4UIcmdWithAnInteger("/rdecay01/histo/printHisto",this);
  prhistoCmd->SetGuidance("print histo #id on ascii file");
  prhistoCmd->SetParameterName("id",false);
  prhistoCmd->SetRange("id>0");
    
  rmhistoCmd = new G4UIcmdWithAnInteger("/rdecay01/histo/removeHisto",this);
  rmhistoCmd->SetGuidance("desactivate histo  #id");
  rmhistoCmd->SetParameterName("id",false);
  rmhistoCmd->SetRange("id>0");
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

HistoMessenger::~HistoMessenger()
{
  delete rmhistoCmd;
  delete prhistoCmd;  
  delete histoCmd;
  delete optionCmd;
  delete typeCmd;  
  delete factoryCmd;
  delete histoDir;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void HistoMessenger::SetNewValue(G4UIcommand* command, G4String newValues)
{
  if (command == factoryCmd)
    histoManager->SetFileName(newValues);

  if (command == typeCmd)
    histoManager->SetFileType(newValues);
    
  if (command == optionCmd)
    histoManager->SetFileOption(newValues);

  if (command == histoCmd)
   { G4int ih,nbBins; G4double vmin,vmax; char unts[30];
     const char* t = newValues;
     std::istringstream is(t);
     is >> ih >> nbBins >> vmin >> vmax >> unts;
     G4String unit = unts;
     G4double vUnit = 1. ;
     if (unit != "none") vUnit = G4UIcommand::ValueOf(unit);
     histoManager->SetHisto (ih,nbBins,vmin*vUnit,vmax*vUnit,unit);
   }
   
  if (command == prhistoCmd)
    histoManager->PrintHisto(prhistoCmd->GetNewIntValue(newValues));
        
  if (command == rmhistoCmd)
    histoManager->RemoveHisto(rmhistoCmd->GetNewIntValue(newValues));         
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

