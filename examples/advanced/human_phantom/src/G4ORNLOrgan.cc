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
// Authors: S. Guatelli and M. G. Pia, INFN Genova, Italy
// 
// Based on code developed by the undergraduate student G. Guerrieri 
// Note: this is a preliminary beta-version of the code; an improved 
// version will be distributed in the next Geant4 public release, compliant
// with the design in a forthcoming publication, and subject to a 
// design and code review.
//
#include "G4ORNLOrgan.hh"
#include "G4PhysicalVolumeStore.hh"
#include "G4Processor/GDMLProcessor.h"
#include "globals.hh"
#include "G4SDManager.hh"
#include "G4VisAttributes.hh"
#include "G4HumanPhantomColour.hh"

G4ORNLOrgan::G4ORNLOrgan()
{
}

G4ORNLOrgan::~G4ORNLOrgan()
{
  sxp.Finalize();
}

G4VPhysicalVolume* G4ORNLOrgan::ConstructOrgan(G4VPhysicalVolume* mother, G4bool sensitivity, G4String gdmlFile, 
G4String colourName, G4bool wireFrame)
{
  G4cout <<"G4ORNLOrgan::ConstructOrgan" << G4endl;
  // Initialize GDML Processor
  sxp.Initialize();
  config.SetURI("gdmlData/" + gdmlFile + ".gdml");
  config.SetSetupName( "Default" );
  sxp.Configure( &config );

  // Run GDML Processor
  sxp.Run(); 
  
  G4int stringLenght = gdmlFile.size();
  G4int i = gdmlFile.find("Female",0);
  //G4cout <<"i: "<< i << G4endl;

  G4int j = gdmlFile.find("Male",0);
  //G4cout <<"j: "<< j << G4endl;

  G4int stringPosition=0;

  if (j == -1) stringPosition = 11;
  else if (i == -1) stringPosition = 9;
  G4String name = gdmlFile.substr(stringPosition, stringLenght - stringPosition);
  //G4cout << name<< G4endl;
  G4String logicalVolumeName = name + "Volume"; 
  G4LogicalVolume* logicOrgan = (G4LogicalVolume *)GDMLProcessor::GetInstance()->GetLogicalVolume(logicalVolumeName);

  G4ThreeVector position = (G4ThreeVector)*GDMLProcessor::GetInstance()->GetPosition("OrganPos");
  G4RotationMatrix* rm = (G4RotationMatrix*)GDMLProcessor::GetInstance()->GetRotation("OrganRot");
 
  G4PhysicalVolumeStore::DeRegister((G4VPhysicalVolume*)GDMLProcessor::GetInstance()->GetWorldVolume());
  // Define rotation and position here!
  G4VPhysicalVolume* physOrgan = new G4PVPlacement(rm,position,
      			       "physicalOrgan",
  			       logicOrgan,
			       mother,
			       false,
			       0);

  // Sensitive Body Part
  if (sensitivity==true)
  { 
    G4SDManager* SDman = G4SDManager::GetSDMpointer();
    logicOrgan->SetSensitiveDetector( SDman->FindSensitiveDetector("BodyPartSD") );
  }
// Visualization Attributes
  G4HumanPhantomColour* colourPointer = new G4HumanPhantomColour();
  G4Colour colour = colourPointer -> GetColour(colourName);
  G4VisAttributes* organVisAtt = new G4VisAttributes(colour);
  // Visualization Attributes
  organVisAtt->SetForceSolid(wireFrame);
  logicOrgan->SetVisAttributes(organVisAtt);

  G4cout << "Organ created !!!!!!  from " <<gdmlFile <<G4endl;
  
  return physOrgan;
}
