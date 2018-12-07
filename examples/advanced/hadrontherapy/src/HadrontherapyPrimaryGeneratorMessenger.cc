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
// Hadrontherapy advanced example for Geant4
// See more at: https://twiki.cern.ch/twiki/bin/view/Geant4/AdvancedExamplesHadrontherapy

#include "HadrontherapyPrimaryGeneratorMessenger.hh"
#include "HadrontherapyPrimaryGeneratorAction.hh"
#include "G4UIdirectory.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"
#include "G4UIcmdWithADouble.hh"
#include "G4SystemOfUnits.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIcmdWithAnInteger.hh"
#include "G4UIcmdWithABool.hh"

HadrontherapyPrimaryGeneratorMessenger::HadrontherapyPrimaryGeneratorMessenger(
                                                                               HadrontherapyPrimaryGeneratorAction* HadrontherapyGun)
:HadrontherapyAction(HadrontherapyGun)
{
    
    changeTheSource = new G4UIdirectory("/primaryParticleData/");
    changeTheSource -> SetGuidance("Command to change the source");
    
    NewSource = new G4UIcmdWithABool("/primaryParticleData/NewSource",this);
    NewSource -> SetParameterName("NewSource ", true);
    NewSource -> SetDefaultValue("false");
    NewSource -> SetGuidance("Set if you want read a .txt file"
                             "\n[usage]: /primaryParticleData/NewSource [true/false]");
    NewSource -> AvailableForStates(G4State_Idle, G4State_PreInit);
    
    
    calculatedPhaseSpaceFileIN=new G4UIcmdWithAString("/primaryParticleData/calculatedPhaseSpaceFileIN",this);
    calculatedPhaseSpaceFileIN->SetDefaultValue("");
    calculatedPhaseSpaceFileIN->SetGuidance("full path and file name of the phase space file to be used as particle generator");
    
    
    
}

HadrontherapyPrimaryGeneratorMessenger::~HadrontherapyPrimaryGeneratorMessenger()
{
    
    delete NewSource;
    delete changeTheSource;
    delete calculatedPhaseSpaceFileIN;
    
    
}


void HadrontherapyPrimaryGeneratorMessenger::SetNewValue(G4UIcommand* command, G4String NewValue)

{
    if (command == calculatedPhaseSpaceFileIN)
    {
        HadrontherapyAction->setCalculatedPhaseSpaceFileIN(NewValue);
    }
    
    if (command == NewSource)
    {
        HadrontherapyAction->setNewSource(NewSource->GetNewBoolValue(NewValue));
        
    }
    
}



