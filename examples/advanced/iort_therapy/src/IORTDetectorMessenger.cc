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
// This is the *BASIC* version of IORT, a Geant4-based application
//
// Main Authors: G.Russo(a,b), C.Casarino*(c), G.C. Candiano(c), G.A.P. Cirrone(d), F.Romano(d)
// Contributor Authors: S.Guatelli(e)
// Past Authors: G.Arnetta(c), S.E.Mazzaglia(d)
//    
//   (a) Fondazione Istituto San Raffaele G.Giglio, Cefalù, Italy
//   (b) IBFM-CNR , Segrate (Milano), Italy
//   (c) LATO (Laboratorio di Tecnologie Oncologiche), Cefalù, Italy
//   (d) Laboratori Nazionali del Sud of the INFN, Catania, Italy
//   (e) University of Wollongong, Australia
//
//   *Corresponding author, email to carlo.casarino@polooncologicocefalu.it
//////////////////////////////////////////////////////////////////////////////////////////////

#include "IORTDetectorMessenger.hh"
#include "IORTDetectorConstruction.hh"
#include "G4UIdirectory.hh"
#include "G4UIcmdWith3VectorAndUnit.hh"
#include "G4UIcmdWithoutParameter.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"   // aggiunto

/////////////////////////////////////////////////////////////////////////////
IORTDetectorMessenger::IORTDetectorMessenger(IORTDetectorConstruction* detector)
  :iortDetector(detector)
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


    // Change disc1 
    changeTheDisc1Dir = new G4UIdirectory("/ProtectionDisc1/");
    changeTheDisc1Dir -> SetGuidance("Command to change the Disc1");
    
    changeOuterRadiusDiscoIORTCmd = new G4UIcmdWithADoubleAndUnit("/ProtectionDisc1/OuterRadiusDisc1",this);
    changeOuterRadiusDiscoIORTCmd -> SetGuidance("Set size of outer radius");
    changeOuterRadiusDiscoIORTCmd -> SetParameterName("Size",false);
    changeOuterRadiusDiscoIORTCmd -> SetDefaultUnit("mm");  
    changeOuterRadiusDiscoIORTCmd -> SetUnitCandidates("mm cm m");  
    changeOuterRadiusDiscoIORTCmd -> AvailableForStates(G4State_Idle);

    changeinnerRadiusDiscoIORTCmd = new G4UIcmdWithADoubleAndUnit("/ProtectionDisc1/InnerRadiusDisc1",this);
    changeinnerRadiusDiscoIORTCmd -> SetGuidance("Set size of inner radius");
    changeinnerRadiusDiscoIORTCmd -> SetParameterName("Size",false);
    changeinnerRadiusDiscoIORTCmd -> SetDefaultUnit("mm");  
    changeinnerRadiusDiscoIORTCmd -> SetUnitCandidates("mm cm m");  
    changeinnerRadiusDiscoIORTCmd -> AvailableForStates(G4State_Idle);


    changeheightDiscoIORTCmd = new G4UIcmdWithADoubleAndUnit("/ProtectionDisc1/HeightDisc1",this);
    changeheightDiscoIORTCmd -> SetGuidance("Set size of higth");
    changeheightDiscoIORTCmd -> SetParameterName("Size",false);
    changeheightDiscoIORTCmd -> SetDefaultUnit("mm");  
    changeheightDiscoIORTCmd -> SetUnitCandidates("mm cm m");  
    changeheightDiscoIORTCmd -> AvailableForStates(G4State_Idle);

    changeDiscoXPositionIORTCmd = new G4UIcmdWithADoubleAndUnit("/ProtectionDisc1/XPositionDisc1",this);
    changeDiscoXPositionIORTCmd -> SetGuidance("Set the X position");
    changeDiscoXPositionIORTCmd -> SetParameterName("Size",false);
    changeDiscoXPositionIORTCmd -> SetDefaultUnit("mm");  
    changeDiscoXPositionIORTCmd -> SetUnitCandidates("mm cm m");  
    changeDiscoXPositionIORTCmd -> AvailableForStates(G4State_Idle);

    changeDiscoYPositionIORTCmd = new G4UIcmdWithADoubleAndUnit("/ProtectionDisc1/YPositionDisc1",this);
    changeDiscoYPositionIORTCmd -> SetGuidance("Set the Y position");
    changeDiscoYPositionIORTCmd -> SetParameterName("Size",false);
    changeDiscoYPositionIORTCmd -> SetDefaultUnit("mm");  
    changeDiscoYPositionIORTCmd -> SetUnitCandidates("mm cm m");  
    changeDiscoYPositionIORTCmd -> AvailableForStates(G4State_Idle);

    changeDiscoZPositionIORTCmd = new G4UIcmdWithADoubleAndUnit("/ProtectionDisc1/ZPositionDisc1",this);
    changeDiscoZPositionIORTCmd -> SetGuidance("Set the Z position");
    changeDiscoZPositionIORTCmd -> SetParameterName("Size",false);
    changeDiscoZPositionIORTCmd -> SetDefaultUnit("mm");  
    changeDiscoZPositionIORTCmd -> SetUnitCandidates("mm cm m");  
    changeDiscoZPositionIORTCmd -> AvailableForStates(G4State_Idle);

    changeTheDisc1MaterialCmd = new G4UIcmdWithAString("/ProtectionDisc1/material", this);
    changeTheDisc1MaterialCmd -> SetGuidance("Change the Disc1 material"); 
    changeTheDisc1MaterialCmd -> SetParameterName("Disc1Material", false);
    changeTheDisc1MaterialCmd -> SetDefaultValue("G4_WATER");
    changeTheDisc1MaterialCmd -> AvailableForStates(G4State_Idle);



    // Change disc2 
    changeTheDisc2Dir = new G4UIdirectory("/ProtectionDisc2/");
    changeTheDisc2Dir -> SetGuidance("Command to change the Disc2");
    
    changeOuterRadiusDisco1IORTCmd = new G4UIcmdWithADoubleAndUnit("/ProtectionDisc2/OuterRadiusDisc2",this);
    changeOuterRadiusDisco1IORTCmd -> SetGuidance("Set size of outer radius");
    changeOuterRadiusDisco1IORTCmd -> SetParameterName("Size",false);
    changeOuterRadiusDisco1IORTCmd -> SetDefaultUnit("mm");  
    changeOuterRadiusDisco1IORTCmd -> SetUnitCandidates("mm cm m");  
    changeOuterRadiusDisco1IORTCmd -> AvailableForStates(G4State_Idle);

    changeinnerRadiusDisco1IORTCmd = new G4UIcmdWithADoubleAndUnit("/ProtectionDisc2/InnerRadiusDisc2",this);
    changeinnerRadiusDisco1IORTCmd -> SetGuidance("Set size of inner radius");
    changeinnerRadiusDisco1IORTCmd -> SetParameterName("Size",false);
    changeinnerRadiusDisco1IORTCmd -> SetDefaultUnit("mm");  
    changeinnerRadiusDisco1IORTCmd -> SetUnitCandidates("mm cm m");  
    changeinnerRadiusDisco1IORTCmd -> AvailableForStates(G4State_Idle);


    changeheightDisco1IORTCmd = new G4UIcmdWithADoubleAndUnit("/ProtectionDisc2/HeightDisc2",this);
    changeheightDisco1IORTCmd -> SetGuidance("Set size of higth");
    changeheightDisco1IORTCmd -> SetParameterName("Size",false);
    changeheightDisco1IORTCmd -> SetDefaultUnit("mm");  
    changeheightDisco1IORTCmd -> SetUnitCandidates("mm cm m");  
    changeheightDisco1IORTCmd -> AvailableForStates(G4State_Idle);

    changeDisco1XPositionIORTCmd = new G4UIcmdWithADoubleAndUnit("/ProtectionDisc2/XPositionDisc2",this);
    changeDisco1XPositionIORTCmd -> SetGuidance("Set the X position");
    changeDisco1XPositionIORTCmd -> SetParameterName("Size",false);
    changeDisco1XPositionIORTCmd -> SetDefaultUnit("mm");  
    changeDisco1XPositionIORTCmd -> SetUnitCandidates("mm cm m");  
    changeDisco1XPositionIORTCmd -> AvailableForStates(G4State_Idle);

    changeTheDisc2MaterialCmd = new G4UIcmdWithAString("/ProtectionDisc2/material", this);
    changeTheDisc2MaterialCmd -> SetGuidance("Change the Disc2 material"); 
    changeTheDisc2MaterialCmd -> SetParameterName("Disc1Material", false);
    changeTheDisc2MaterialCmd -> SetDefaultValue("G4_WATER");
    changeTheDisc2MaterialCmd -> AvailableForStates(G4State_Idle);

    // Delete disc 1-2
    deleteTheDiscDir = new G4UIdirectory("/DeleteProtectionDisc/");
    deleteTheDiscDir -> SetGuidance("Command to delete the 1-2 Discs ");

    deletediscCmd = new G4UIcmdWithoutParameter("/DeleteProtectionDisc/delete",this);
    deletediscCmd->SetGuidance("Delete the 1-2 Discs geometry.");
    deletediscCmd->SetGuidance("This command MUST be applied before \"beamOn\" ");
    deletediscCmd->AvailableForStates(G4State_Idle);

    // Insert disc 1-2
    insertTheDiscDir = new G4UIdirectory("/InsertProtectionDisc/");
    insertTheDiscDir -> SetGuidance("Command to insert the 1-2 Discs ");

    insertdiscCmd = new G4UIcmdWithoutParameter("/InsertProtectionDisc/insert",this);
    insertdiscCmd->SetGuidance("Insert the 1-2 Discs geometry.");
    insertdiscCmd->SetGuidance("This command MUST be applied before \"beamOn\" ");
    insertdiscCmd->SetGuidance("After this command MUST be applied update command \"beamOn\" ");
    insertdiscCmd->AvailableForStates(G4State_Idle);

    // Change Tilt angle disc1 + disc2
    changeTheAnglediscDir = new G4UIdirectory("/ChangeTiltAngleDisc1-2/");
    changeTheAnglediscDir -> SetGuidance("Set tilt angle of the 1-2 Discs");

    changeTheAnglediscCmd = new G4UIcmdWithADoubleAndUnit("/ChangeTiltAngleDisc1-2/TiltAngleDisc1-2",this);
    changeTheAnglediscCmd -> SetParameterName("Angle",false);
    changeTheAnglediscCmd -> SetDefaultUnit("deg");  
    changeTheAnglediscCmd -> SetUnitCandidates("deg rad");  
    changeTheAnglediscCmd -> AvailableForStates(G4State_Idle); 
   }

/////////////////////////////////////////////////////////////////////////////
IORTDetectorMessenger::~IORTDetectorMessenger()
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
    
    delete changeTheDisc1Dir;
    delete changeOuterRadiusDiscoIORTCmd;
    delete changeinnerRadiusDiscoIORTCmd;
    delete changeheightDiscoIORTCmd; 
    delete changeDiscoXPositionIORTCmd;
    delete changeDiscoYPositionIORTCmd;
    delete changeDiscoZPositionIORTCmd;
    delete changeTheDisc1MaterialCmd;

    delete changeTheDisc2Dir;
    delete changeOuterRadiusDisco1IORTCmd;
    delete changeinnerRadiusDisco1IORTCmd;
    delete changeheightDisco1IORTCmd; 
    delete changeDisco1XPositionIORTCmd;
    delete changeTheDisc2MaterialCmd;

    delete deletediscCmd;
    delete insertdiscCmd;

    delete changeTheAnglediscDir;
    delete changeTheAnglediscCmd;
}

/////////////////////////////////////////////////////////////////////////////
void IORTDetectorMessenger::SetNewValue(G4UIcommand* command,G4String newValue)
{
	
  if( command == changeThePhantomSizeCmd)
    {
	G4ThreeVector size = changeThePhantomSizeCmd -> GetNew3VectorValue(newValue);
	iortDetector -> SetPhantomSize(size.getX(),size.getY(),size.getZ());
    }
  else if (command == changeThePhantomPositionCmd )
  {
	 G4ThreeVector size = changeThePhantomPositionCmd -> GetNew3VectorValue(newValue);
         iortDetector -> SetPhantomPosition(size);
  }
  else if (command == changeThePhantomMaterialCmd)
  {
      iortDetector -> SetPhantomMaterial(newValue);
  }
  else if (command == changeTheDetectorSizeCmd)
  {
	G4ThreeVector size = changeTheDetectorSizeCmd  -> GetNew3VectorValue(newValue);
        iortDetector -> SetDetectorSize(size.getX(),size.getY(),size.getZ());
  }
  else if (command == changeTheDetectorToPhantomPositionCmd)
  {
	G4ThreeVector size = changeTheDetectorToPhantomPositionCmd-> GetNew3VectorValue(newValue);
        iortDetector -> SetDetectorToPhantomPosition(size);
  }
  else if (command == changeTheDetectorVoxelCmd)
  {
	G4ThreeVector size = changeTheDetectorVoxelCmd  -> GetNew3VectorValue(newValue);
        iortDetector -> SetVoxelSize(size.getX(),size.getY(),size.getZ());
  }
/////////////////////disc/////////////////////////////////

  else if( command == changeOuterRadiusDiscoIORTCmd )
    { 
         iortDetector -> SetOuterRadiusDiscoIORT
	(changeOuterRadiusDiscoIORTCmd -> GetNewDoubleValue(newValue));
    }
  else if( command == changeinnerRadiusDiscoIORTCmd )
    { iortDetector -> SetinnerRadiusDiscoIORT
	(changeinnerRadiusDiscoIORTCmd -> GetNewDoubleValue(newValue));
    }
  else if( command == changeheightDiscoIORTCmd )
    { iortDetector -> SetheightDiscoIORT
	(changeheightDiscoIORTCmd -> GetNewDoubleValue(newValue));
    }
  else if( command == changeDiscoXPositionIORTCmd )
    { iortDetector -> SetDiscoXPositionIORT
	(changeDiscoXPositionIORTCmd -> GetNewDoubleValue(newValue));
    }
  else if( command == changeDiscoYPositionIORTCmd )
    { iortDetector -> SetDiscoYPositionIORT
	(changeDiscoYPositionIORTCmd -> GetNewDoubleValue(newValue));
    }
  else if( command == changeDiscoZPositionIORTCmd )
    { iortDetector -> SetDiscoZPositionIORT
	(changeDiscoZPositionIORTCmd -> GetNewDoubleValue(newValue));
    }
  else if (command == changeTheDisc1MaterialCmd)
  {
      iortDetector -> SetDiscoMaterialIORT(newValue);
  }

  else if( command == changeOuterRadiusDisco1IORTCmd )
    { iortDetector -> SetOuterRadiusDiscoIORT1
	(changeOuterRadiusDisco1IORTCmd -> GetNewDoubleValue(newValue));
    }
  else if( command == changeinnerRadiusDisco1IORTCmd )
    { iortDetector -> SetinnerRadiusDiscoIORT1
	(changeinnerRadiusDisco1IORTCmd -> GetNewDoubleValue(newValue));
    }
  else if( command == changeheightDisco1IORTCmd )
    { iortDetector -> SetheightDiscoIORT1
	(changeheightDisco1IORTCmd -> GetNewDoubleValue(newValue));
    }
  else if( command == changeDisco1XPositionIORTCmd )
    { iortDetector -> SetDiscoXPositionIORT1
	(changeDisco1XPositionIORTCmd -> GetNewDoubleValue(newValue));
    }
  else if (command == changeTheDisc2MaterialCmd)
  {
      iortDetector -> SetDiscoMaterialIORT1(newValue);
  }
  else if (command == changeTheAnglediscCmd)
  {
      iortDetector -> SetAngleDiscoIORT0
        (changeTheAnglediscCmd -> GetNewDoubleValue(newValue));
  }  



  else if (command == updateCmd)
  {
      iortDetector -> UpdateGeometry();
  }

  else if (command == deletediscCmd)
  {
      iortDetector -> DeleteDisc();
  }

else if (command == insertdiscCmd)
  {
      iortDetector -> ConstructDisc();
  }
}
