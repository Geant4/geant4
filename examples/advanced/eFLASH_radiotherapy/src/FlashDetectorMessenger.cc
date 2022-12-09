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
//
//////////////////////////////////////////////////////////////////////////////////////////////

#include "FlashDetectorMessenger.hh"
#include "FlashDetectorConstruction.hh"
#include "G4UIdirectory.hh"
#include "G4UIcmdWith3VectorAndUnit.hh"
#include "G4UIcmdWithoutParameter.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"  

/////////////////////////////////////////////////////////////////////////////
FlashDetectorMessenger::FlashDetectorMessenger(FlashDetectorConstruction* detector)
  :flashDetector(detector)
{
    // fChange Phantom size
    fChangeThePhantomDir = new G4UIdirectory("/fChangePhantom/");
    fChangeThePhantomDir -> SetGuidance("Command to fChange the Phantom Size/position");

     // fChange Phantom material 
    fChangeThePhantomMaterialCmd = new G4UIcmdWithAString("/changePhantom/material", this);
    fChangeThePhantomMaterialCmd -> SetGuidance("fChange the Phantom and the detector material"); 
    fChangeThePhantomMaterialCmd -> SetParameterName("PhantomMaterial", false);
    fChangeThePhantomMaterialCmd -> SetDefaultValue("G4_WATER");
    fChangeThePhantomMaterialCmd -> AvailableForStates(G4State_Idle);

   
 
    // fChange Detector size
    fChangeTheDetectorDir = new G4UIdirectory("/changeDetector/");
    fChangeTheDetectorDir -> SetGuidance("Command to fChange the Detector Size");

// fChange Detector material 
    fChangeTheDetectorMaterialCmd = new G4UIcmdWithAString("/changeDetector/material", this);
    fChangeTheDetectorMaterialCmd  -> SetGuidance("fChange the Phantom and the detector material"); 
    fChangeTheDetectorMaterialCmd  -> SetParameterName("PhantomMaterial", false);
   
    fChangeTheDetectorMaterialCmd  -> AvailableForStates(G4State_Idle);

    

    fUpdateCmd = new G4UIcmdWithoutParameter("/changePhantom/update",this);
    fUpdateCmd->SetGuidance("Update Phantom/Detector geometry.");
    fUpdateCmd->SetGuidance("This command MUST be applied before \"beamOn\" ");
    fUpdateCmd->SetGuidance("if you fChanged geometrical value(s).");
    fUpdateCmd->AvailableForStates(G4State_Idle);

    fUpdateCmd_d = new G4UIcmdWithoutParameter("/changeDetector/update",this);
    fUpdateCmd_d->SetGuidance("Update Detector geometry.");
    fUpdateCmd_d->SetGuidance("This command MUST be applied before \"beamOn\" ");
    fUpdateCmd_d->SetGuidance("if you fChanged geometrical value(s).");
    fUpdateCmd_d->AvailableForStates(G4State_Idle);
    


   }

/////////////////////////////////////////////////////////////////////////////
FlashDetectorMessenger::~FlashDetectorMessenger()
{


  delete fChangeThePhantomDir;
    delete fChangeTheDetectorDir;
  
    delete fChangeThePhantomMaterialCmd;
    delete fChangeTheDetectorMaterialCmd; 
    
}

/////////////////////////////////////////////////////////////////////////////
void FlashDetectorMessenger::SetNewValue(G4UIcommand* command,G4String newValue)
{
	
 
  if (command == fChangeThePhantomMaterialCmd)
  {
      flashDetector -> SetPhantomMaterial(newValue);
  }

    else if (command == fChangeTheDetectorMaterialCmd)
  {
      flashDetector -> SetDetectorMaterial(newValue);
  }
  
  
  
}
