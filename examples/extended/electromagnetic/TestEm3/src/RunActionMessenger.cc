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
// $Id: RunActionMessenger.cc,v 1.6 2004-06-09 14:18:48 maire Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "RunActionMessenger.hh"

#include "RunAction.hh"
#include "G4UIcommand.hh"
#include "G4UIparameter.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIcmdWith3Vector.hh"

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

  accCmd1 = new G4UIcmdWith3Vector("/testem/run/acceptanceL1",this);
  accCmd1->SetGuidance("set Edep and RMS");
  accCmd1->SetGuidance("acceptance values for first layer");
  accCmd1->SetParameterName("edep","rms","limit",true);
  accCmd1->SetRange("edep>0 && edep<1 && rms>0");
  accCmd1->AvailableForStates(G4State_PreInit,G4State_Idle);

  accCmd2 = new G4UIcmdWith3Vector("/testem/run/acceptanceL2",this);
  accCmd2->SetGuidance("set Edep and RMS");
  accCmd2->SetGuidance("acceptance values for 2nd layer");
  accCmd2->SetParameterName("edep","rms","limit",true);
  accCmd2->SetRange("edep>0 && edep<1 && rms>0");
  accCmd2->AvailableForStates(G4State_PreInit,G4State_Idle);  
    
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

RunActionMessenger::~RunActionMessenger()
{
  delete HistoCmd;  
  delete fileCmd;
  delete accCmd1;
  delete accCmd2;  
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
   
  if( command == accCmd1 )
   { Run->SetEdepAndRMS(0,accCmd1->GetNew3VectorValue(newValue));}

  if( command == accCmd2 )
   { Run->SetEdepAndRMS(1,accCmd2->GetNew3VectorValue(newValue));}
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
