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
#include "G4ORNLLeftKidney.hh"
#include "G4PhysicalVolumeStore.hh"
#include "G4Processor/GDMLProcessor.h"
#include "globals.hh"
#include "G4SDManager.hh"

#include "G4VisAttributes.hh"

G4ORNLLeftKidney::G4ORNLLeftKidney()
{
}

G4ORNLLeftKidney::~G4ORNLLeftKidney()
{
  sxp.Finalize();
}

G4VPhysicalVolume* G4ORNLLeftKidney::ConstructLeftKidney(G4VPhysicalVolume* mother, G4String sex, G4bool sensitivity)
{
  // Initialize GDML Processor
  sxp.Initialize();
  config.SetURI( "gdmlData/"+sex+"/ORNLLeftKidney.gdml" );
  config.SetSetupName( "Default" );
  sxp.Configure( &config );

  // Run GDML Processor
  sxp.Run();
 
  G4cout<< "Left kidney, gdml file read!!" << G4endl;

  G4LogicalVolume* logicLeftKidney = (G4LogicalVolume *)GDMLProcessor::GetInstance()->GetLogicalVolume("LeftKidneyVolume");

  G4ThreeVector position = (G4ThreeVector)*GDMLProcessor::GetInstance()->GetPosition("LeftKidneyPos");
  G4RotationMatrix* rm = (G4RotationMatrix*)GDMLProcessor::GetInstance()->GetRotation("LeftKidneyRot");
  G4cout << "positiona nd rotation read!" << G4endl;
G4PhysicalVolumeStore::DeRegister((G4VPhysicalVolume*)GDMLProcessor::GetInstance()->GetWorldVolume());  
  // Define rotation and position here!
  G4VPhysicalVolume* physLeftKidney = new G4PVPlacement(rm,position,
      			       "physicalLeftKidney",
  			       logicLeftKidney,
			       mother,
			       false,
			       0);


  // Sensitive Body Part
  if (sensitivity==true)
  { 
    G4SDManager* SDman = G4SDManager::GetSDMpointer();
    logicLeftKidney->SetSensitiveDetector( SDman->FindSensitiveDetector("BodyPartSD") );
  }

  // Visualization Attributes
  G4VisAttributes* LeftKidneyVisAtt = new G4VisAttributes(G4Colour(0.72,0.52,0.04));
  LeftKidneyVisAtt->SetForceSolid(true);
  logicLeftKidney->SetVisAttributes(LeftKidneyVisAtt);

  G4cout << "LeftKidney created !!!!!!" << G4endl;
  
  return physLeftKidney;
}
