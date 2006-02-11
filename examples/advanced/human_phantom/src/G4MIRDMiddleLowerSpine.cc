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
#include "G4MIRDMiddleLowerSpine.hh"

#include "G4Processor/GDMLProcessor.h"
#include "globals.hh"
#include "G4SDManager.hh"
#include "G4VisAttributes.hh"

G4MIRDMiddleLowerSpine::G4MIRDMiddleLowerSpine()
{
}

G4MIRDMiddleLowerSpine::~G4MIRDMiddleLowerSpine()
{
  sxp.Finalize();
}

G4VPhysicalVolume* G4MIRDMiddleLowerSpine::ConstructMiddleLowerSpine(G4VPhysicalVolume* mother, G4String sex, G4bool sensitivity)
{
  // Initialize GDML Processor
  sxp.Initialize();
  config.SetURI( "gdmlData/"+sex+"/MIRDMiddleLowerSpine.gdml" );
  config.SetSetupName( "Default" );
  sxp.Configure( &config );

  // Run GDML Processor
  sxp.Run();
 

  G4LogicalVolume* logicMiddleLowerSpine = (G4LogicalVolume *)GDMLProcessor::GetInstance()->GetLogicalVolume("MiddleLowerSpineVolume");

  G4ThreeVector position = (G4ThreeVector)*GDMLProcessor::GetInstance()->GetPosition("MiddleLowerSpinePos");
  G4RotationMatrix* rm = (G4RotationMatrix*)GDMLProcessor::GetInstance()->GetRotation("MiddleLowerSpineRot");
  
  // Define rotation and position here!
  G4VPhysicalVolume* physMiddleLowerSpine = new G4PVPlacement(rm,position,
      			       "physicalMiddleLowerSpine",
  			       logicMiddleLowerSpine,
			       mother,
			       false,
			       0);

  // Sensitive Body Part
  if (sensitivity==true)
  { 
    G4SDManager* SDman = G4SDManager::GetSDMpointer();
    logicMiddleLowerSpine->SetSensitiveDetector( SDman->FindSensitiveDetector("BodyPartSD") );
  }

  // Visualization Attributes
  G4VisAttributes* MiddleLowerSpineVisAtt = new G4VisAttributes(G4Colour(0.46,0.53,0.6));
  MiddleLowerSpineVisAtt->SetForceSolid(true);
  logicMiddleLowerSpine->SetVisAttributes(MiddleLowerSpineVisAtt);

  G4cout << "MiddleLowerSpine created !!!!!!" << G4endl;

  // Testing MiddleLowerSpine Volume
  G4double MiddleLowerSpineVol = logicMiddleLowerSpine->GetSolid()->GetCubicVolume();
  G4cout << "Volume of MiddleLowerSpine = " << MiddleLowerSpineVol/cm3 << " cm^3" << G4endl;
  
  // Testing MiddleLowerSpine Material
  G4String MiddleLowerSpineMat = logicMiddleLowerSpine->GetMaterial()->GetName();
  G4cout << "Material of MiddleLowerSpine = " << MiddleLowerSpineMat << G4endl;
  
  // Testing Density
  G4double MiddleLowerSpineDensity = logicMiddleLowerSpine->GetMaterial()->GetDensity();
  G4cout << "Density of Material = " << MiddleLowerSpineDensity*cm3/g << " g/cm^3" << G4endl;

  // Testing Mass
  G4double MiddleLowerSpineMass = (MiddleLowerSpineVol)*MiddleLowerSpineDensity;
  G4cout << "Mass of MiddleLowerSpine = " << MiddleLowerSpineMass/gram << " g" << G4endl;

  
  return physMiddleLowerSpine;
}
