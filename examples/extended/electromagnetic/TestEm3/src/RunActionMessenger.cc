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
// $Id: RunActionMessenger.cc,v 1.3 2004/01/21 17:29:27 maire Exp $
// GEANT4 tag $Name: geant4-06-00-patch-01 $
//
// 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "RunActionMessenger.hh"

#include "RunAction.hh"
#include "G4UIcommand.hh"
#include "G4UIparameter.hh"
#include "G4UIcmdWithAString.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

RunActionMessenger::RunActionMessenger(RunAction* run)
:Run(run)
{    
  fileCmd = new G4UIcmdWithAString("/testem/run/setFileName",this);
  fileCmd->SetGuidance("set name for the histograms file");

  HistoCmd = new G4UIcommand("/testem/run/setHisto",this);
  HistoCmd->SetGuidance("Set histo Edep in absorber k");
  HistoCmd->SetGuidance("  histo=absor number : from 0 to NbOfAbsor-1");
  HistoCmd->SetGuidance("  nb bins; Emin; Emax; unit (of Emin and Emax");
  //
  G4UIparameter* AbsNbPrm = new G4UIparameter("AbsorNb",'i',false);
  AbsNbPrm->SetGuidance("histo=absor number : from 0 to NbOfAbsor-1");
  AbsNbPrm->SetParameterRange("AbsorNb>=0");
  HistoCmd->SetParameter(AbsNbPrm);
  //
  G4UIparameter* NbinPrm = new G4UIparameter("Nbin",'i',false);
  NbinPrm->SetGuidance("number of bins");
  NbinPrm->SetParameterRange("Nbin>=0");
  HistoCmd->SetParameter(NbinPrm);  
  //    
  G4UIparameter* VminPrm = new G4UIparameter("Vmin",'d',false);
  VminPrm->SetGuidance("Emin, expressed in unit");
  HistoCmd->SetParameter(VminPrm);
  //    
  G4UIparameter* VmaxPrm = new G4UIparameter("Vmax",'d',false);
  VmaxPrm->SetGuidance("Emax, expressed in unit");
  HistoCmd->SetParameter(VmaxPrm);
  //
  G4UIparameter* UnitPrm = new G4UIparameter("unit",'s',false);
  HistoCmd->SetParameter(UnitPrm);  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

RunActionMessenger::~RunActionMessenger()
{
  delete HistoCmd;  
  delete fileCmd;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void RunActionMessenger::SetNewValue(G4UIcommand* command,G4String newValue)
{   
  if (command == fileCmd) Run->SetFileName(newValue);

  if (command == HistoCmd)
   { G4int idh,nbins; G4double vmin,vmax; char unts[30];
     const char* t = newValue;
     std::istrstream is((char*)t);
     is >> idh >> nbins >> vmin >> vmax >> unts;
     G4String unit  = unts;
     G4double vUnit = G4UIcommand::ValueOf(unit);
     Run->SetHisto (idh,nbins,vmin*vUnit,vmax*vUnit,unit);
   }         
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
