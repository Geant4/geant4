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
#include "G4MIRDSkull.hh"

#include "G4Processor/GDMLProcessor.h"
#include "globals.hh"
#include "G4SDManager.hh"
#include "G4VisAttributes.hh"

G4MIRDSkull::G4MIRDSkull()
{
}

G4MIRDSkull::~G4MIRDSkull()
{
  sxp.Finalize();
}

G4VPhysicalVolume* G4MIRDSkull::ConstructSkull(G4VPhysicalVolume* mother, G4String sex, G4bool sensitivity)
{
  // Initialize GDML Processor
  sxp.Initialize();
  config.SetURI( "gdmlData/"+sex+"/MIRDSkull.gdml" );
  config.SetSetupName( "Default" );
  sxp.Configure( &config );

  // Run GDML Processor
  sxp.Run();
 

  G4LogicalVolume* logicSkull = (G4LogicalVolume *)GDMLProcessor::GetInstance()->GetLogicalVolume("SkullVolume");

  G4ThreeVector position = (G4ThreeVector)*GDMLProcessor::GetInstance()->GetPosition("SkullPos");
  G4RotationMatrix* rm = (G4RotationMatrix*)GDMLProcessor::GetInstance()->GetRotation("SkullRot");
  
  // Define rotation and position here!
  G4VPhysicalVolume* physSkull = new G4PVPlacement(rm,position,
      			       "physicalSkull",
  			       logicSkull,
			       mother,
			       false,
			       0);

  // Sensitive Body Part
  if (sensitivity==true)
  { 
    G4SDManager* SDman = G4SDManager::GetSDMpointer();
    logicSkull->SetSensitiveDetector( SDman->FindSensitiveDetector("BodyPartSD") );
  }

  // Visualization Attributes
  G4VisAttributes* SkullVisAtt = new G4VisAttributes(G4Colour(0.46,0.53,0.6));
  SkullVisAtt->SetForceSolid(true);
  logicSkull->SetVisAttributes(SkullVisAtt);

  G4cout << "Skull created !!!!!!" << G4endl;


  // Testing Skull Volume
  G4double SkullVol = logicSkull->GetSolid()->GetCubicVolume();
  G4cout << "Volume of Skull = " << SkullVol/cm3 << " cm^3" << G4endl;
  
  // Testing Skull Material
  G4String SkullMat = logicSkull->GetMaterial()->GetName();
  G4cout << "Material of Skull = " << SkullMat << G4endl;
  
  // Testing Density
  G4double SkullDensity = logicSkull->GetMaterial()->GetDensity();
  G4cout << "Density of Material = " << SkullDensity*cm3/g << " g/cm^3" << G4endl;

  // Testing Mass
  G4double SkullMass = (SkullVol)*SkullDensity;
  G4cout << "Mass of Skull = " << SkullMass/gram << " g" << G4endl;

  
  return physSkull;
}
