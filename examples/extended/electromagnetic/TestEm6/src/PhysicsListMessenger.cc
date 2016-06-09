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
// $Id: PhysicsListMessenger.cc,v 1.8 2006/06/29 16:57:11 gunter Exp $
// GEANT4 tag $Name: geant4-09-01 $
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
 
  GammaToMuPairFac=new G4UIcmdWithADouble
                                      ("/testem/phys/SetGammaToMuPairFac",this);
  GammaToMuPairFac->SetGuidance(
         "Set factor to artificially increase the GammaToMuPair cross section");
  GammaToMuPairFac->SetParameterName("GammaToMuPairFac",false);
  GammaToMuPairFac->SetRange("GammaToMuPairFac>0.0");
  GammaToMuPairFac->AvailableForStates(G4State_PreInit,G4State_Idle);

  AnnihiToMuPairFac=new G4UIcmdWithADouble
                                     ("/testem/phys/SetAnnihiToMuPairFac",this);
  AnnihiToMuPairFac->SetGuidance(
        "Set factor to artificially increase the AnnihiToMuPair cross section");
  AnnihiToMuPairFac->SetParameterName("AnnihiToMuPairFac",false);
  AnnihiToMuPairFac->SetRange("AnnihiToMuPairFac>0.0");
  AnnihiToMuPairFac->AvailableForStates(G4State_PreInit,G4State_Idle);
  
  AnnihiToHadronFac=new G4UIcmdWithADouble
                                      ("/testem/phys/SetAnnihiToHadronFac",this);
  AnnihiToHadronFac->SetGuidance(
       "Set factor to artificially increase the AnnihiToHadrons cross section");
  AnnihiToHadronFac->SetParameterName("AnnihiToHadFac",false);
  AnnihiToHadronFac->SetRange("AnnihiToHadFac>0.0");
  AnnihiToHadronFac->AvailableForStates(G4State_PreInit,G4State_Idle);  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PhysicsListMessenger::~PhysicsListMessenger()
{
  delete GammaToMuPairFac;
  delete AnnihiToMuPairFac;
  delete AnnihiToHadronFac;  
  delete physDir;  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PhysicsListMessenger::SetNewValue(G4UIcommand* command,
                                          G4String newValue)
{ 
  if(command == GammaToMuPairFac)
   { physList->SetGammaToMuPairFac(
                           GammaToMuPairFac->GetNewDoubleValue(newValue));}   

  if( command == AnnihiToMuPairFac)
   { physList->SetAnnihiToMuPairFac(
                          AnnihiToMuPairFac->GetNewDoubleValue(newValue));}
			  
  if( command == AnnihiToHadronFac)
   { physList->SetAnnihiToHadronFac(
                          AnnihiToHadronFac->GetNewDoubleValue(newValue));}			     
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
