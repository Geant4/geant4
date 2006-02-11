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
#include "G4MIRDArmBone.hh"

#include "G4Processor/GDMLProcessor.h"
#include "globals.hh"
#include "G4SDManager.hh"
#include "G4VisAttributes.hh"

G4MIRDArmBone::G4MIRDArmBone()
{
}

G4MIRDArmBone::~G4MIRDArmBone()
{
  sxp.Finalize();
}

G4VPhysicalVolume* G4MIRDArmBone::ConstructArmBone(G4VPhysicalVolume* mother, G4String sex, G4bool sensitivity)
{
  // Initialize GDML Processor
  sxp.Initialize();
  config.SetURI( "gdmlData/"+sex+"/MIRDArmBone.gdml" );
  config.SetSetupName( "Default" );
  sxp.Configure( &config );

  // Run GDML Processor
  sxp.Run();
 

  G4LogicalVolume* logicArmBone = (G4LogicalVolume *)GDMLProcessor::GetInstance()->GetLogicalVolume("ArmBoneVolume");

  G4ThreeVector position = (G4ThreeVector)*GDMLProcessor::GetInstance()->GetPosition("ArmBonePos");
  G4RotationMatrix* rm = (G4RotationMatrix*)GDMLProcessor::GetInstance()->GetRotation("ArmBoneRot");
  
  // Define rotation and position here!
  G4VPhysicalVolume* physArmBone = new G4PVPlacement(rm,position,
      			       "physicalArmBone",
  			       logicArmBone,
			       mother,
			       false,
			       0);

  // Sensitive Body Part
  if (sensitivity==true)
  { 
    G4SDManager* SDman = G4SDManager::GetSDMpointer();
    logicArmBone->SetSensitiveDetector( SDman->FindSensitiveDetector("BodyPartSD") );
  }


  // Visualization Attributes
  G4VisAttributes* ArmBoneVisAtt = new G4VisAttributes(G4Colour(0.46,0.53,0.6));
  ArmBoneVisAtt->SetForceSolid(true);
  logicArmBone->SetVisAttributes(ArmBoneVisAtt);

  G4cout << "ArmBone created !!!!!!" << G4endl;
 
  // Testing ArmBone Volume
  G4double ArmBoneVol = logicArmBone->GetSolid()->GetCubicVolume();
  G4cout << "Volume of ArmBone = " << ArmBoneVol/cm3 << " cm^3" << G4endl;
  
  // Testing ArmBone Material
  G4String ArmBoneMat = logicArmBone->GetMaterial()->GetName();
  G4cout << "Material of ArmBone = " << ArmBoneMat << G4endl;
  
  // Testing Density
  G4double ArmBoneDensity = logicArmBone->GetMaterial()->GetDensity();
  G4cout << "Density of Material = " << ArmBoneDensity*cm3/g << " g/cm^3" << G4endl;

  // Testing Mass
  G4double ArmBoneMass = (ArmBoneVol)*ArmBoneDensity;
  G4cout << "Mass of ArmBone = " << ArmBoneMass/gram << " g" << G4endl;

  
  return physArmBone;
}
