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
// $Id: Em3RunActionMessenger.cc,v 1.7 2001-10-22 10:58:59 maire Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "Em3RunActionMessenger.hh"

#include "Em3RunAction.hh"
#include "G4UIdirectory.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIcmdWithAnInteger.hh"
#include "G4UIcommand.hh"
#include "G4UIparameter.hh"
#include "G4ios.hh"
#include "globals.hh"
#include "Randomize.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

Em3RunActionMessenger::Em3RunActionMessenger(Em3RunAction* run)
:Em3Run(run)
{    
  RndmDir = new G4UIdirectory("/rndm/");
  RndmDir->SetGuidance("Rndm status control.");
  
  RndmSaveCmd = new G4UIcmdWithAnInteger("/rndm/save",this);
  RndmSaveCmd->SetGuidance("set frequency to save rndm on external files.");
  RndmSaveCmd->SetGuidance("freq = 0 not saved");
  RndmSaveCmd->SetGuidance("freq > 0 saved on: beginOfRun.rndm");
  RndmSaveCmd->SetGuidance("freq > 0 saved on:   endOfRun.rndm");
  RndmSaveCmd->SetGuidance("freq = 2 saved on: beginOfEvent.rndm");    
  RndmSaveCmd->SetParameterName("frequency",false);
  RndmSaveCmd->SetRange("frequency>=0 && frequency<=2");
  RndmSaveCmd->AvailableForStates(PreInit,Idle); 
         
  RndmReadCmd = new G4UIcmdWithAString("/rndm/read",this);
  RndmReadCmd->SetGuidance("get rndm status from an external file.");
  RndmReadCmd->SetParameterName("fileName",true);
  RndmReadCmd->SetDefaultValue ("beginOfRun.rndm");
  RndmReadCmd->AvailableForStates(PreInit,Idle);
  
  HistoCmd = new G4UIcommand("/run/setHisto",this);
  HistoCmd->SetGuidance("Set histo Edep/Ebeam in absorber k");
  HistoCmd->SetGuidance("  histo=absor number : from 0 to NbOfAbsor-1");
  HistoCmd->SetGuidance("  number of bins; Emin/Ebeam; Emax/Ebeam");
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
  VminPrm->SetGuidance("Vmin=Emin/Ebeam");
  VminPrm->SetParameterRange("Vmin>=0.&&Vmin<=1.");
  HistoCmd->SetParameter(VminPrm);
  //    
  G4UIparameter* VmaxPrm = new G4UIparameter("Vmax",'d',false);
  VmaxPrm->SetGuidance("Vmax=Emax/Ebeam");
  VmaxPrm->SetParameterRange("Vmax>=0.&&Vmax<=1.");
  HistoCmd->SetParameter(VmaxPrm);  
  //
  HistoCmd->AvailableForStates(Idle);    
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

Em3RunActionMessenger::~Em3RunActionMessenger()
{
  delete RndmSaveCmd; delete RndmReadCmd; delete RndmDir;
  delete HistoCmd;  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void Em3RunActionMessenger::SetNewValue(G4UIcommand* command,G4String newValue)
{   
  if (command == RndmSaveCmd)
      Em3Run->SetRndmFreq(RndmSaveCmd->GetNewIntValue(newValue));
		 
  if (command == RndmReadCmd)
    { G4cout << "\n---> rndm status restored from file: " << newValue << G4endl;
      HepRandom::restoreEngineStatus(newValue);
      HepRandom::showEngineStatus();
    }
    
  if (command == HistoCmd)
   { G4int idh,nbins; G4double vmin,vmax;
     const char* t = newValue;
     G4std::istrstream is((char*)t);
     is >> idh >> nbins >> vmin >> vmax;
     Em3Run->SetHisto (idh,nbins,vmin,vmax);
   }         
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
