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
// $Id: PhysicsListMessenger.cc,v 1.4 2003-01-27 09:46:20 hbu Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "PhysicsListMessenger.hh"

#include "PhysicsList.hh"
#include "G4UIcmdWithADouble.hh"
#include "globals.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PhysicsListMessenger::PhysicsListMessenger(PhysicsList* physL)
:physList(physL)
{ 
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
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PhysicsListMessenger::~PhysicsListMessenger()
{
  delete GammaToMuPairFac;
  delete AnnihiToMuPairFac;
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
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
