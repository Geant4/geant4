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
/// \file NeutronHPMessenger.cc
/// \brief Implementation of the NeutronHPMessenger class
//
// $Id: NeutronHPMessenger.cc 67268 2013-02-13 11:38:40Z ihrivnac $
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "NeutronHPMessenger.hh"

#include "NeutronHPphysics.hh"

#include "G4UIdirectory.hh"
#include "G4UIcmdWithABool.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

NeutronHPMessenger::NeutronHPMessenger(NeutronHPphysics* phys)
:G4UImessenger(),fNeutronPhysics(phys),
 fPhysDir(0), fThermalCmd(0)
{ 
  fPhysDir = new G4UIdirectory("/testhadr/phys/");
  fPhysDir->SetGuidance("physics list commands");
   
  fThermalCmd = new G4UIcmdWithABool("/testhadr/phys/thermalScattering",this);
  fThermalCmd->SetGuidance("set thermal scattering model");
  fThermalCmd->SetParameterName("thermal",false);
  fThermalCmd->AvailableForStates(G4State_PreInit);  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

NeutronHPMessenger::~NeutronHPMessenger()
{
  delete fThermalCmd;
  delete fPhysDir;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void NeutronHPMessenger::SetNewValue(G4UIcommand* command, G4String newValue)
{   
  if (command == fThermalCmd)
   {fNeutronPhysics->SetThermalPhysics(fThermalCmd->GetNewBoolValue(newValue));}
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
