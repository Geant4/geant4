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
// $Id: EmAnalysisMessenger.cc,v 1.1 2004-08-26 15:07:50 vnivanch Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "EmAnalysisMessenger.hh"

#include "EmAnalysis.hh"
#include "G4UIdirectory.hh"
#include "G4UIcommand.hh"
#include "G4UIparameter.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIcmdWithAnInteger.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

EmAnalysisMessenger::EmAnalysisMessenger(EmAnalysis* manager)
:ema (manager)
{
  emaDir = new G4UIdirectory("/testem/calc/");
  emaDir->SetGuidance("cross sections calculation");

  saveCmd = new G4UIcommand("/testem/calc/save",this);
  saveCmd->SetGuidance("Fill histograms and save to file");

  verbCmd = new G4UIcmdWithAnInteger("/testem/calc/verbose",this);
  verbCmd->SetGuidance("set verbose level");

  listCmd = new G4UIcmdWithAnInteger("/testem/calc/list",this);
  listCmd->SetGuidance("printout of histogram parameters");

  histoCmd = new G4UIcommand("/testem/calc/hist",this);
  histoCmd->SetGuidance("Set new histogram:");
  histoCmd->SetGuidance(" particle; material; process; type; histId");
  //
  G4UIparameter* particle = new G4UIparameter("particle",'s',false);
  particle->SetGuidance("particle name");
  histoCmd->SetParameter(particle);
  //
  G4UIparameter* material = new G4UIparameter("material",'s',false);
  material->SetGuidance("material name");
  histoCmd->SetParameter(particle);
  //
  G4UIparameter* process = new G4UIparameter("process",'s',false);
  process->SetGuidance("process name");
  histoCmd->SetParameter(process);
  //
  G4UIparameter* type = new G4UIparameter("type",'s',false);
  type->SetGuidance("type: dedx or cross");
  type->SetDefaultValue("dedx");
  histoCmd->SetParameter(type);
  //
  G4UIparameter* id = new G4UIparameter("histId",'s',false);
  id->SetGuidance("histogram id");
  id->SetDefaultValue("");
  histoCmd->SetParameter(id);
  //

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

EmAnalysisMessenger::~EmAnalysisMessenger()
{
  delete histoCmd;
  delete verbCmd;
  delete listCmd;
  //  delete setCmd;
  //  delete actCmd;
  delete saveCmd;
  delete emaDir;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void EmAnalysisMessenger::SetNewValue(G4UIcommand* command,G4String newValues)
{
  if (command == saveCmd)
    ema->saveToFile();

  if (command == verbCmd)
    ema->setVerbose(verbCmd->GetNewIntValue(newValues));

  if (command == listCmd)
    ema->PrintHist(listCmd->GetNewIntValue(newValues));
    
  if (command == histoCmd) { 
    G4String part, mat, proc, type, hid;
    const char* t = newValues;
    std::istrstream is((char*)t);
    is >> part >> mat >> proc >> type >> hid;
    G4int n = ema->AddHistOnCrossSection(part,mat,proc,type,hid);
    G4cout << "### Histogram ID= " << n << " is booked" << G4endl;
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
