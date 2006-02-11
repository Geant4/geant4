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
#include "G4MIRDLung.hh"

#include "G4Processor/GDMLProcessor.h"
#include "globals.hh"
#include "G4SDManager.hh"
#include "G4VisAttributes.hh"

G4MIRDLung::G4MIRDLung()
{
}

G4MIRDLung::~G4MIRDLung()
{
  sxp.Finalize();
}

G4VPhysicalVolume* G4MIRDLung::ConstructLung(G4VPhysicalVolume* mother, G4String sex, G4bool sensitivity)
{
  // Initialize GDML Processor
  sxp.Initialize();
  config.SetURI( "gdmlData/"+sex+"/MIRDLung.gdml" );
  config.SetSetupName( "Default" );
  sxp.Configure( &config );

  // Run GDML Processor
  sxp.Run();
 

  G4LogicalVolume* logicLung = (G4LogicalVolume *)GDMLProcessor::GetInstance()->GetLogicalVolume("LungVolume");

  G4ThreeVector position = (G4ThreeVector)*GDMLProcessor::GetInstance()->GetPosition("LungPos");
  G4RotationMatrix* rm = (G4RotationMatrix*)GDMLProcessor::GetInstance()->GetRotation("LungRot");
  
  // Define rotation and position here!
  G4VPhysicalVolume* physLung = new G4PVPlacement(rm,position,
      			       "physicalLung",
  			       logicLung,
			       mother,
			       false,
			       0);

  // Sensitive Body Part
  if (sensitivity==true)
  { 
    G4SDManager* SDman = G4SDManager::GetSDMpointer();
    logicLung->SetSensitiveDetector( SDman->FindSensitiveDetector("BodyPartSD") );
  }

  // Visualization Attributes
  G4VisAttributes* LungVisAtt = new G4VisAttributes(G4Colour(0.25,0.41,0.88));
  LungVisAtt->SetForceSolid(true);
  logicLung->SetVisAttributes(LungVisAtt);

  G4cout << "Lung created !!!!!!" << G4endl;

  // Testing Lung Volume
  G4double LungVol = logicLung->GetSolid()->GetCubicVolume();
  G4cout << "Volume of Lung = " << LungVol/cm3 << " cm^3" << G4endl;
  
  // Testing Lung Material
  G4String LungMat = logicLung->GetMaterial()->GetName();
  G4cout << "Material of Lung = " << LungMat << G4endl;
  
  // Testing Density
  G4double LungDensity = logicLung->GetMaterial()->GetDensity();
  G4cout << "Density of Material = " << LungDensity*cm3/g << " g/cm^3" << G4endl;

  // Testing Mass
  G4double LungMass = (LungVol)*LungDensity;
  G4cout << "Mass of Lung = " << LungMass/gram << " g" << G4endl;
  
  return physLung;
}
