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
//
// $Id: HistoMessenger.cc,v 1.1 2003/08/11 10:20:44 maire Exp $
// GEANT4 tag $Name: geant4-06-00 $
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#ifdef G4ANALYSIS_USE

#include "HistoMessenger.hh"

#include "HistoManager.hh"
#include "G4UIdirectory.hh"
#include "G4UIcommand.hh"
#include "G4UIparameter.hh"
#include "G4UIcmdWithAString.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

HistoMessenger::HistoMessenger(HistoManager* manager)
:histoManager (manager)
{
  histoDir = new G4UIdirectory("/testem/histo/");
  histoDir->SetGuidance("histograms control");

  factoryCmd = new G4UIcmdWithAString("/testem/histo/setFileName",this);
  factoryCmd->SetGuidance("set name for the histograms file");
   
  histoCmd = new G4UIcommand("/testem/histo/setHisto",this);
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
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

HistoMessenger::~HistoMessenger()
{
  delete histoCmd;
  delete factoryCmd;
  delete histoDir;  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void HistoMessenger::SetNewValue(G4UIcommand* command,G4String newValues)
{
  if (command == factoryCmd)
    histoManager->SetFactory(newValues);
    
  if (command == histoCmd)
   { G4int ih,nbBins; G4double vmin,vmax; char unts[30];
     const char* t = newValues;
     std::istrstream is((char*)t);
     is >> ih >> nbBins >> vmin >> vmax >> unts;
     G4String unit = unts;
     G4double vUnit = 1. ;
     if (unit != "none") vUnit = G4UIcommand::ValueOf(unit);
     histoManager->SetHisto (ih,nbBins,vmin*vUnit,vmax*vUnit,unit);
   }      
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif      
