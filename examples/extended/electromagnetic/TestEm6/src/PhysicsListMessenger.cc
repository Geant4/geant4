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
// $Id: PhysicsListMessenger.cc,v 1.9 2009-11-27 14:54:58 hbu Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "PhysicsListMessenger.hh"

#include "PhysicsList.hh"
#include "G4UIdirectory.hh"
#include "G4UIcmdWithADouble.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PhysicsListMessenger::PhysicsListMessenger(PhysicsList* physL)
:physList(physL)
{
  physDir = new G4UIdirectory("/testem/phys/");
  physDir->SetGuidance("physics list commands");
 
  GammaToMuPairFacCmd=new G4UIcmdWithADouble
                                      ("/testem/phys/SetGammaToMuPairFac",this);
  GammaToMuPairFacCmd->SetGuidance(
         "Set factor to artificially increase the GammaToMuPair cross section");
  GammaToMuPairFacCmd->SetParameterName("GammaToMuPairFac",false);
  GammaToMuPairFacCmd->SetRange("GammaToMuPairFac>0.0");
  GammaToMuPairFacCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  AnnihiToMuPairFacCmd=new G4UIcmdWithADouble
                                     ("/testem/phys/SetAnnihiToMuPairFac",this);
  AnnihiToMuPairFacCmd->SetGuidance(
        "Set factor to artificially increase the AnnihiToMuPair cross section");
  AnnihiToMuPairFacCmd->SetParameterName("AnnihiToMuPairFac",false);
  AnnihiToMuPairFacCmd->SetRange("AnnihiToMuPairFac>0.0");
  AnnihiToMuPairFacCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
  
  AnnihiToHadronFacCmd=new G4UIcmdWithADouble
                                      ("/testem/phys/SetAnnihiToHadronFac",this);
  AnnihiToHadronFacCmd->SetGuidance(
       "Set factor to artificially increase the AnnihiToHadrons cross section");
  AnnihiToHadronFacCmd->SetParameterName("AnnihiToHadFac",false);
  AnnihiToHadronFacCmd->SetRange("AnnihiToHadFac>0.0");
  AnnihiToHadronFacCmd->AvailableForStates(G4State_PreInit,G4State_Idle);  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PhysicsListMessenger::~PhysicsListMessenger()
{
  delete GammaToMuPairFacCmd;
  delete AnnihiToMuPairFacCmd;
  delete AnnihiToHadronFacCmd;  
  delete physDir;  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PhysicsListMessenger::SetNewValue(G4UIcommand* command,
                                          G4String newValue)
{ 
  if(command == GammaToMuPairFacCmd)
   { physList->SetGammaToMuPairFac(
                           GammaToMuPairFacCmd->GetNewDoubleValue(newValue));}   

  if( command == AnnihiToMuPairFacCmd)
   { physList->SetAnnihiToMuPairFac(
                          AnnihiToMuPairFacCmd->GetNewDoubleValue(newValue));}
			  
  if( command == AnnihiToHadronFacCmd)
   { physList->SetAnnihiToHadronFac(
                          AnnihiToHadronFacCmd->GetNewDoubleValue(newValue));}			     
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
