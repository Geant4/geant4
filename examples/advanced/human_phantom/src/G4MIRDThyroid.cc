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
#include "G4MIRDThyroid.hh"

#include "G4Processor/GDMLProcessor.h"
#include "globals.hh"
#include "G4SDManager.hh"
#include "G4VisAttributes.hh"

G4MIRDThyroid::G4MIRDThyroid()
{
}

G4MIRDThyroid::~G4MIRDThyroid()
{
  sxp.Finalize();
}

G4VPhysicalVolume* G4MIRDThyroid::ConstructThyroid(G4VPhysicalVolume* mother, G4String sex, G4bool sensitivity)
{
  // Initialize GDML Processor
  sxp.Initialize();
  config.SetURI( "gdmlData/"+sex+"/MIRDThyroid.gdml" );
  config.SetSetupName( "Default" );
  sxp.Configure( &config );

  // Run GDML Processor
  sxp.Run();
 

  G4LogicalVolume* logicThyroid = (G4LogicalVolume *)GDMLProcessor::GetInstance()->GetLogicalVolume("ThyroidVolume");

  G4ThreeVector position = (G4ThreeVector)*GDMLProcessor::GetInstance()->GetPosition("ThyroidPos");
  G4RotationMatrix* rm = (G4RotationMatrix*)GDMLProcessor::GetInstance()->GetRotation("ThyroidRot");
  
  // Define rotation and position here!
  G4VPhysicalVolume* physThyroid = new G4PVPlacement(rm,position,
      			       "physicalThyroid",
  			       logicThyroid,
			       mother,
			       false,
			       0);

  // Sensitive Body Part
  if (sensitivity==true)
  { 
    G4SDManager* SDman = G4SDManager::GetSDMpointer();
    logicThyroid->SetSensitiveDetector( SDman->FindSensitiveDetector("BodyPartSD") );
  }

  // Visualization Attributes
  G4VisAttributes* ThyroidVisAtt = new G4VisAttributes(G4Colour(0.0,1.0,1.0));
  ThyroidVisAtt->SetForceSolid(true);
  logicThyroid->SetVisAttributes(ThyroidVisAtt);

  G4cout << "Thyroid created !!!!!!" << G4endl;

  // Testing Thyroid Volume
  G4double ThyroidVol = logicThyroid->GetSolid()->GetCubicVolume();
  G4cout << "Volume of Thyroid = " << ThyroidVol/cm3 << " cm^3" << G4endl;
  
  // Testing Thyroid Material
  G4String ThyroidMat = logicThyroid->GetMaterial()->GetName();
  G4cout << "Material of Thyroid = " << ThyroidMat << G4endl;
  
  // Testing Density
  G4double ThyroidDensity = logicThyroid->GetMaterial()->GetDensity();
  G4cout << "Density of Material = " << ThyroidDensity*cm3/g << " g/cm^3" << G4endl;

  // Testing Mass
  G4double ThyroidMass = (ThyroidVol)*ThyroidDensity;
  G4cout << "Mass of Thyroid = " << ThyroidMass/gram << " g" << G4endl;

  
  return physThyroid;
}
