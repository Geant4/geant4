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
#include "G4MIRDBreast.hh"

#include "G4Processor/GDMLProcessor.h"
#include "globals.hh"
#include "G4SDManager.hh"
#include "G4VisAttributes.hh"

G4MIRDBreast::G4MIRDBreast()
{
}

G4MIRDBreast::~G4MIRDBreast()
{
  sxp.Finalize();
}

G4VPhysicalVolume* G4MIRDBreast::ConstructBreast(G4VPhysicalVolume* mother, G4String sex, G4bool sensitivity)
{
  // Initialize GDML Processor
  sxp.Initialize();
  config.SetURI( "gdmlData/"+sex+"/MIRDBreast.gdml" );
  config.SetSetupName( "Default" );
  sxp.Configure( &config );

  // Run GDML Processor
  sxp.Run();
 

  G4LogicalVolume* logicBreast = (G4LogicalVolume *)GDMLProcessor::GetInstance()->GetLogicalVolume("BreastVolume");

  G4ThreeVector position = (G4ThreeVector)*GDMLProcessor::GetInstance()->GetPosition("BreastPos");
  G4RotationMatrix* rm = (G4RotationMatrix*)GDMLProcessor::GetInstance()->GetRotation("BreastRot");
  
  // Define rotation and position here!
  G4VPhysicalVolume* physBreast = new G4PVPlacement(rm,position,
      			       "physicalBreast",
  			       logicBreast,
			       mother,
			       false,
			       0);

  // Sensitive Body Part
  if (sensitivity==true)
  { 
    G4SDManager* SDman = G4SDManager::GetSDMpointer();
    logicBreast->SetSensitiveDetector( SDman->FindSensitiveDetector("BodyPartSD") );
  }

  // Visualization Attributes
  G4VisAttributes* BreastVisAtt = new G4VisAttributes(G4Colour(1.0,0.41,0.71));
  BreastVisAtt->SetForceSolid(true);
  logicBreast->SetVisAttributes(BreastVisAtt);

  G4cout << "Breast created !!!!!!" << G4endl;
  
  // Testing Breast Volume
  G4double BreastVol = logicBreast->GetSolid()->GetCubicVolume();
  G4cout << "Volume of Breast = " << BreastVol/cm3 << " cm^3" << G4endl;
  
  // Testing Breast Material
  G4String BreastMat = logicBreast->GetMaterial()->GetName();
  G4cout << "Material of Breast = " << BreastMat << G4endl;
  
  // Testing Density
  G4double BreastDensity = logicBreast->GetMaterial()->GetDensity();
  G4cout << "Density of Material = " << BreastDensity*cm3/g << " g/cm^3" << G4endl;

  // Testing Mass
  G4double BreastMass = (BreastVol)*BreastDensity;
  G4cout << "Mass of Breast = " << BreastMass/gram << " g" << G4endl;


  return physBreast;
}
