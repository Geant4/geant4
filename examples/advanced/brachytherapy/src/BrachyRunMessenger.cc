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
// $Id: BrachyRunMessenger.cc,v 1.1 2002-11-15 17:31:23 guatelli Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "BrachyRunMessenger.hh"
#include"BrachyFactoryIr.hh"
#include "BrachyRunAction.hh"
#include "G4UIdirectory.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIcmdWithAnInteger.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"
#include "G4UIcmdWithoutParameter.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

BrachyRunMessenger::BrachyRunMessenger( BrachyRunAction* pBrachyRun):
 pRun(pBrachyRun)
{ 

  mydetDir = new G4UIdirectory("/run/");
  mydetDir->SetGuidance("Control gamma energy.");
  selDetCmd=new G4UIcmdWithAString("/run/energy",this);
  
 
  
  selDetCmd->SetGuidance("Select the energy of gamma emitted by the source.");
  selDetCmd->SetGuidance(" Iodium:  Iodium Source ");
  selDetCmd->SetGuidance("  Iridium: Iridium  Source  "  );
  selDetCmd->SetParameterName("choice",true);
  selDetCmd->SetDefaultValue("Iridium");
  selDetCmd->SetCandidates("Iridium / Iodium");
  selDetCmd->AvailableForStates(PreInit,Idle);

 
 }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

BrachyRunMessenger::~BrachyRunMessenger()
{
 
  delete selDetCmd; 
  delete   mydetDir;
  //delete detDir;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void BrachyRunMessenger::SetNewValue(G4UIcommand* command,G4String newValue)
{ 
  

 if( command == selDetCmd )
  {if (newValue== "Iodium") pRun->SelectEnergy(1);
    else  pRun->SelectEnergy(0);
  }
   
}

