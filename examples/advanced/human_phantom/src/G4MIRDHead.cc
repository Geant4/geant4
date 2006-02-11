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
#include "G4MIRDHead.hh"

#include "G4Processor/GDMLProcessor.h"
#include "globals.hh"

#include "G4SDManager.hh"

#include "G4VisAttributes.hh"

G4MIRDHead::G4MIRDHead()
{
}

G4MIRDHead::~G4MIRDHead()
{
  sxp.Finalize();
}

G4VPhysicalVolume* G4MIRDHead::ConstructHead(G4VPhysicalVolume* mother, G4String sex, G4bool sensitivity)
{
  // Initialize GDML Processor
  sxp.Initialize();
  config.SetURI( "gdmlData/"+sex+"/MIRDHead.gdml" );
  config.SetSetupName( "Default" );
  sxp.Configure( &config );

  // Run GDML Processor
  sxp.Run();
 

  G4LogicalVolume* logicHead = (G4LogicalVolume *)GDMLProcessor::GetInstance()->GetLogicalVolume("HeadVolume");

  G4ThreeVector position = (G4ThreeVector)*GDMLProcessor::GetInstance()->GetPosition("HeadPos");
  G4RotationMatrix* rm = (G4RotationMatrix*)GDMLProcessor::GetInstance()->GetRotation("HeadRot");
  
  // Define rotation and position here!
  G4VPhysicalVolume* physHead = new G4PVPlacement(rm,position,
      			       "physicalHead",
  			       logicHead,
			       mother,
			       false,
			       0);

  // Sensitive Body Part
  if (sensitivity==true)
  { 
    G4SDManager* SDman = G4SDManager::GetSDMpointer();
    logicHead->SetSensitiveDetector( SDman->FindSensitiveDetector("BodyPartSD") );
  }

  // Visualization Attributes
  G4VisAttributes* HeadVisAtt = new G4VisAttributes(G4Colour(0.94,0.5,0.5));
  HeadVisAtt->SetForceSolid(false);
  logicHead->SetVisAttributes(HeadVisAtt);

  G4cout << "Head created !!!!!!" << G4endl;

  // Testing Head Volume
  G4double HeadVol = logicHead->GetSolid()->GetCubicVolume();
  G4cout << "Volume of Head = " << HeadVol/cm3 << " cm^3" << G4endl;
  
  // Testing Head Material
  G4String HeadMat = logicHead->GetMaterial()->GetName();
  G4cout << "Material of Head = " << HeadMat << G4endl;
  
  // Testing Density
  G4double HeadDensity = logicHead->GetMaterial()->GetDensity();
  G4cout << "Density of Material = " << HeadDensity*cm3/g << " g/cm^3" << G4endl;

  // Testing Mass
  G4double HeadMass = (HeadVol)*HeadDensity;
  G4cout << "Mass of Head = " << HeadMass/gram << " g" << G4endl;


  
  return physHead;
}
