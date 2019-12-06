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
//  Author: F. Poignant, floriane.poignant@gmail.com
//

#include "STCyclotronPrimaryGeneratorActionMessenger.hh"
#include "STCyclotronPrimaryGeneratorAction.hh"
#include "G4UIdirectory.hh"
#include "G4UIcmdWithADouble.hh"
#include "G4SystemOfUnits.hh"
#include "G4UnitsTable.hh"

STCyclotronPrimaryGeneratorActionMessenger::STCyclotronPrimaryGeneratorActionMessenger(STCyclotronPrimaryGeneratorAction* primary)
  :fG4Primary(primary)
{
    // Change beam current
    fBeamCurrent = new G4UIdirectory("/setBeamCurrent/");
    fBeamCurrent -> SetGuidance("Change the beam current of the cyclotron");

    fChangeBeamCurrentCmd = new G4UIcmdWithADouble("/setBeamCurrent/beamCurrent", this);
    fChangeBeamCurrentCmd -> SetGuidance("Change the value of the current (in ampere)."
	                                   "\nThe default value is 30E-6 ampere.");
    fChangeBeamCurrentCmd -> SetParameterName("BeamCurrent", true);
    fChangeBeamCurrentCmd -> SetRange("BeamCurrent > 0.");
    fChangeBeamCurrentCmd -> SetDefaultValue(30.E-6);
    fChangeBeamCurrentCmd -> AvailableForStates(G4State_Idle);
    
   }

STCyclotronPrimaryGeneratorActionMessenger::~STCyclotronPrimaryGeneratorActionMessenger()
{
    delete fBeamCurrent; 
    delete fChangeBeamCurrentCmd; 
    
}

void STCyclotronPrimaryGeneratorActionMessenger::SetNewValue(G4UIcommand* command,G4String newValue)
{
  
  if( command == fChangeBeamCurrentCmd)
    {
      G4double updatedValue  = fChangeBeamCurrentCmd -> GetNewDoubleValue(newValue);
      fG4Primary -> SetBeamCurrent(updatedValue);
    }
  
}
