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
#include "G4MIRDSpleen.hh"

#include "G4Processor/GDMLProcessor.h"
#include "globals.hh"
#include "G4SDManager.hh"
#include "G4VisAttributes.hh"

G4MIRDSpleen::G4MIRDSpleen()
{
}

G4MIRDSpleen::~G4MIRDSpleen()
{
  sxp.Finalize();
}

G4VPhysicalVolume* G4MIRDSpleen::ConstructSpleen(G4VPhysicalVolume* mother, G4String sex, G4bool sensitivity)
{
  // Initialize GDML Processor
  sxp.Initialize();
  config.SetURI( "gdmlData/"+sex+"/MIRDSpleen.gdml" );
  config.SetSetupName( "Default" );
  sxp.Configure( &config );

  // Run GDML Processor
  sxp.Run();
 

  G4LogicalVolume* logicSpleen = (G4LogicalVolume *)GDMLProcessor::GetInstance()->GetLogicalVolume("SpleenVolume");

  G4ThreeVector position = (G4ThreeVector)*GDMLProcessor::GetInstance()->GetPosition("SpleenPos");
  G4RotationMatrix* rm = (G4RotationMatrix*)GDMLProcessor::GetInstance()->GetRotation("SpleenRot");
  
  // Define rotation and position here!
  G4VPhysicalVolume* physSpleen = new G4PVPlacement(rm,position,
      			       "physicalSpleen",
  			       logicSpleen,
			       mother,
			       false,
			       0);

  // Sensitive Body Part
  if (sensitivity==true)
  { 
    G4SDManager* SDman = G4SDManager::GetSDMpointer();
    logicSpleen->SetSensitiveDetector( SDman->FindSensitiveDetector("BodyPartSD") );
  }

  // Visualization Attributes
  G4VisAttributes* SpleenVisAtt = new G4VisAttributes(G4Colour(0.41,0.41,0.41));
  SpleenVisAtt->SetForceSolid(true);
  logicSpleen->SetVisAttributes(SpleenVisAtt);

  G4cout << "Spleen created !!!!!!" << G4endl;

  // Testing Spleen Volume
  G4double SpleenVol = logicSpleen->GetSolid()->GetCubicVolume();
  G4cout << "Volume of Spleen = " << SpleenVol/cm3 << " cm^3" << G4endl;
  
  // Testing Spleen Material
  G4String SpleenMat = logicSpleen->GetMaterial()->GetName();
  G4cout << "Material of Spleen = " << SpleenMat << G4endl;
  
  // Testing Density
  G4double SpleenDensity = logicSpleen->GetMaterial()->GetDensity();
  G4cout << "Density of Material = " << SpleenDensity*cm3/g << " g/cm^3" << G4endl;

  // Testing Mass
  G4double SpleenMass = (SpleenVol)*SpleenDensity;
  G4cout << "Mass of Spleen = " << SpleenMass/gram << " g" << G4endl;


  
  return physSpleen;
}
