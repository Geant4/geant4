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
#include "G4MIRDUpperLargeIntestine.hh"

#include "G4Processor/GDMLProcessor.h"
#include "globals.hh"
#include "G4SDManager.hh"
#include "G4VisAttributes.hh"

G4MIRDUpperLargeIntestine::G4MIRDUpperLargeIntestine()
{
}

G4MIRDUpperLargeIntestine::~G4MIRDUpperLargeIntestine()
{
  sxp.Finalize();
}

G4VPhysicalVolume* G4MIRDUpperLargeIntestine::ConstructUpperLargeIntestine(G4VPhysicalVolume* mother, G4String sex, G4bool sensitivity)
{
  // Initialize GDML Processor
  sxp.Initialize();
  config.SetURI( "gdmlData/"+sex+"/MIRDUpperLargeIntestine.gdml" );
  config.SetSetupName( "Default" );
  sxp.Configure( &config );

  // Run GDML Processor
  sxp.Run();
 

  G4LogicalVolume* logicUpperLargeIntestine = (G4LogicalVolume *)GDMLProcessor::GetInstance()->GetLogicalVolume("UpperLargeIntestineVolume");

  G4ThreeVector position = (G4ThreeVector)*GDMLProcessor::GetInstance()->GetPosition("UpperLargeIntestinePos");
  G4RotationMatrix* rm = (G4RotationMatrix*)GDMLProcessor::GetInstance()->GetRotation("UpperLargeIntestineRot");
  
  // Define rotation and position here!
  G4VPhysicalVolume* physUpperLargeIntestine = new G4PVPlacement(rm,position,
      			       "physicalUpperLargeIntestine",
  			       logicUpperLargeIntestine,
			       mother,
			       false,
			       0);

  // Sensitive Body Part
  if (sensitivity==true)
  { 
    G4SDManager* SDman = G4SDManager::GetSDMpointer();
    logicUpperLargeIntestine->SetSensitiveDetector( SDman->FindSensitiveDetector("BodyPartSD") );
  }

  // Visualization Attributes
  G4VisAttributes* UpperLargeIntestineVisAtt = new G4VisAttributes(G4Colour(1.0,1.0,0.0));
  UpperLargeIntestineVisAtt->SetForceSolid(true);
  logicUpperLargeIntestine->SetVisAttributes(UpperLargeIntestineVisAtt);

  G4cout << "UpperLargeIntestine created !!!!!!" << G4endl;

  // Testing UpperLargeIntestine Volume
  G4double UpperLargeIntestineVol = logicUpperLargeIntestine->GetSolid()->GetCubicVolume();
  G4cout << "Volume of UpperLargeIntestine = " << UpperLargeIntestineVol/cm3 << " cm^3" << G4endl;
  
  // Testing UpperLargeIntestine Material
  G4String UpperLargeIntestineMat = logicUpperLargeIntestine->GetMaterial()->GetName();
  G4cout << "Material of UpperLargeIntestine = " << UpperLargeIntestineMat << G4endl;
  
  // Testing Density
  G4double UpperLargeIntestineDensity = logicUpperLargeIntestine->GetMaterial()->GetDensity();
  G4cout << "Density of Material = " << UpperLargeIntestineDensity*cm3/g << " g/cm^3" << G4endl;

  // Testing Mass
  G4double UpperLargeIntestineMass = (UpperLargeIntestineVol)*UpperLargeIntestineDensity;
  G4cout << "Mass of UpperLargeIntestine = " << UpperLargeIntestineMass/gram << " g" << G4endl;

  
  return physUpperLargeIntestine;
}
