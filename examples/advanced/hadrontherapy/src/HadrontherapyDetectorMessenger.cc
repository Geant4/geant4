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

#include "HadrontherapyDetectorMessenger.hh"
#include "HadrontherapyDetectorConstruction.hh"
#include "G4UIdirectory.hh"
#include "G4UIcmdWith3VectorAndUnit.hh"
#include "G4UIcmdWithoutParameter.hh"
#include "G4UIcmdWithAString.hh"
#include "G4SystemOfUnits.hh"
#include "G4UIcmdWithABool.hh"

/////////////////////////////////////////////////////////////////////////////
HadrontherapyDetectorMessenger::HadrontherapyDetectorMessenger(HadrontherapyDetectorConstruction* detector)
  :hadrontherapyDetector(detector)
{
    // Change Phantom size
    changeThePhantomDir = new G4UIdirectory("/changePhantom/");
    changeThePhantomDir -> SetGuidance("Command to change the Phantom Size/position");
    changeThePhantomSizeCmd = new G4UIcmdWith3VectorAndUnit("/changePhantom/size", this);
    changeThePhantomSizeCmd -> SetGuidance("Insert sizes X Y and Z"
	                                   "\n   0 or negative values mean <<Don't change it!>>");
    changeThePhantomSizeCmd -> SetParameterName("PhantomSizeAlongX", 
						"PhantomSizeAlongY", 
						"PhantomSizeAlongZ", false);
    changeThePhantomSizeCmd -> SetDefaultUnit("mm");
    changeThePhantomSizeCmd -> SetUnitCandidates("nm um mm cm"); 
    changeThePhantomSizeCmd -> AvailableForStates(G4State_Idle);


    // Change Phantom material 
    changeThePhantomMaterialCmd = new G4UIcmdWithAString("/changePhantom/material", this);
    changeThePhantomMaterialCmd -> SetGuidance("Change the Phantom and the detector material"); 
    changeThePhantomMaterialCmd -> SetParameterName("PhantomMaterial", false);
    changeThePhantomMaterialCmd -> SetDefaultValue("G4_WATER");
    changeThePhantomMaterialCmd -> AvailableForStates(G4State_Idle);

    // Change Phantom position
    changeThePhantomPositionCmd = new G4UIcmdWith3VectorAndUnit("/changePhantom/position", this);
    changeThePhantomPositionCmd -> SetGuidance("Insert X Y and Z dimensions for the position of the center of the Phantom"
						" respect to that of the \"World\""); 
    changeThePhantomPositionCmd -> SetParameterName("PositionAlongX", 
						    "PositionAlongY", 
						    "PositionAlongZ", false);
    changeThePhantomPositionCmd -> SetDefaultUnit("mm");
    changeThePhantomPositionCmd -> SetUnitCandidates("um mm cm m"); 
    changeThePhantomPositionCmd -> AvailableForStates(G4State_Idle);


    updateCmd = new G4UIcmdWithoutParameter("/changePhantom/update",this);
    updateCmd->SetGuidance("Update Phantom/Detector geometry.");
    updateCmd->SetGuidance("This command MUST be applied before \"beamOn\" ");
    updateCmd->SetGuidance("if you changed geometrical value(s).");
    updateCmd->AvailableForStates(G4State_Idle);

    //  Change detector size
    changeTheDetectorDir = new G4UIdirectory("/changeDetector/");
    changeTheDetectorDir -> SetGuidance("Command to change the Detector's Size/position/Voxels");
    
    changeTheDetectorSizeCmd = new G4UIcmdWith3VectorAndUnit("/changeDetector/size",this);
    changeTheDetectorSizeCmd -> SetGuidance("Insert sizes for X Y and Z dimensions of the Detector"
					    "\n   0 or negative values mean <<Don't change it>>");
    changeTheDetectorSizeCmd -> SetParameterName("DetectorSizeAlongX", "DetectorSizeAlongY", "DetectorSizeAlongZ", false);
    changeTheDetectorSizeCmd -> SetDefaultUnit("mm");
    changeTheDetectorSizeCmd -> SetUnitCandidates("nm um mm cm"); 
    changeTheDetectorSizeCmd -> AvailableForStates(G4State_Idle);

    //  Change the detector to phantom displacement
    changeTheDetectorToPhantomPositionCmd = new G4UIcmdWith3VectorAndUnit("/changeDetector/displacement",this);
    changeTheDetectorToPhantomPositionCmd -> SetGuidance("Insert X Y and Z displacements between Detector and Phantom"
	                                                 "\nNegative values mean <<Don't change it!>>"); 
    changeTheDetectorToPhantomPositionCmd -> SetParameterName("DisplacementAlongX",
							      "DisplacementAlongY", 
							      "DisplacementAlongZ", false);
    changeTheDetectorToPhantomPositionCmd -> SetDefaultUnit("mm");
    changeTheDetectorToPhantomPositionCmd -> SetUnitCandidates("nm um mm cm"); 
    changeTheDetectorToPhantomPositionCmd -> AvailableForStates(G4State_Idle);
    
    // Change voxels by its size
    changeTheDetectorVoxelCmd = new G4UIcmdWith3VectorAndUnit("/changeDetector/voxelSize",this);
    changeTheDetectorVoxelCmd -> SetGuidance("Insert Voxel sizes for X Y and Z dimensions"
		                             "\n   0 or negative values mean <<Don't change it!>>");
    changeTheDetectorVoxelCmd -> SetParameterName("VoxelSizeAlongX", "VoxelSizeAlongY", "VoxelSizeAlongZ", false);
    changeTheDetectorVoxelCmd -> SetDefaultUnit("mm");
    changeTheDetectorVoxelCmd -> SetUnitCandidates("nm um mm cm");
    changeTheDetectorVoxelCmd -> AvailableForStates(G4State_Idle);


    
    

}

/////////////////////////////////////////////////////////////////////////////
HadrontherapyDetectorMessenger::~HadrontherapyDetectorMessenger()
{
    delete changeThePhantomDir; 
    delete changeThePhantomSizeCmd; 
    delete changeThePhantomPositionCmd; 
    delete changeThePhantomMaterialCmd; 
    delete updateCmd;
    delete changeTheDetectorDir; 
    delete changeTheDetectorSizeCmd; 
    delete changeTheDetectorToPhantomPositionCmd; 
    delete changeTheDetectorVoxelCmd;

}

/////////////////////////////////////////////////////////////////////////////
void HadrontherapyDetectorMessenger::SetNewValue(G4UIcommand* command,G4String newValue)
{
	
  if( command == changeThePhantomSizeCmd)
    {
	G4ThreeVector size = changeThePhantomSizeCmd -> GetNew3VectorValue(newValue);
	hadrontherapyDetector -> SetPhantomSize(size.getX(),size.getY(),size.getZ());
    }
  else if (command == changeThePhantomPositionCmd )
  {
	 G4ThreeVector size = changeThePhantomPositionCmd -> GetNew3VectorValue(newValue);
         hadrontherapyDetector -> SetPhantomPosition(size);
  }
  else if (command == changeThePhantomMaterialCmd)
  {
      hadrontherapyDetector -> SetPhantomMaterial(newValue);
  }
  else if (command == changeTheDetectorSizeCmd)
  {
	G4ThreeVector size = changeTheDetectorSizeCmd  -> GetNew3VectorValue(newValue);
        hadrontherapyDetector -> SetDetectorSize(size.getX(),size.getY(),size.getZ());
  }
  else if (command == changeTheDetectorToPhantomPositionCmd)
  {
	G4ThreeVector size = changeTheDetectorToPhantomPositionCmd-> GetNew3VectorValue(newValue);
        hadrontherapyDetector -> SetDetectorToPhantomPosition(size);
  }
  else if (command == changeTheDetectorVoxelCmd)
  {
	G4ThreeVector size = changeTheDetectorVoxelCmd  -> GetNew3VectorValue(newValue);
        hadrontherapyDetector -> SetVoxelSize(size.getX(),size.getY(),size.getZ());
  }
  else if (command == updateCmd)
  {
      hadrontherapyDetector -> UpdateGeometry();
  }
 
}
