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
#include "G4MIRDUrinaryBladder.hh"

#include "G4Processor/GDMLProcessor.h"
#include "globals.hh"
#include "G4SDManager.hh"
#include "G4VisAttributes.hh"

G4MIRDUrinaryBladder::G4MIRDUrinaryBladder()
{
}

G4MIRDUrinaryBladder::~G4MIRDUrinaryBladder()
{
  sxp.Finalize();
}

G4VPhysicalVolume* G4MIRDUrinaryBladder::ConstructUrinaryBladder(G4VPhysicalVolume* mother, G4String sex, G4bool sensitivity)
{
  // Initialize GDML Processor
  sxp.Initialize();
  config.SetURI( "gdmlData/"+sex+"/MIRDUrinaryBladder.gdml" );
  config.SetSetupName( "Default" );
  sxp.Configure( &config );

  // Run GDML Processor
  sxp.Run();
 

  G4LogicalVolume* logicUrinaryBladder = (G4LogicalVolume *)GDMLProcessor::GetInstance()->GetLogicalVolume("UrinaryBladderVolume");

  G4ThreeVector position = (G4ThreeVector)*GDMLProcessor::GetInstance()->GetPosition("UrinaryBladderPos");
  G4RotationMatrix* rm = (G4RotationMatrix*)GDMLProcessor::GetInstance()->GetRotation("UrinaryBladderRot");
  
  // Define rotation and position here!
  G4VPhysicalVolume* physUrinaryBladder = new G4PVPlacement(rm,position,
      			       "physicalUrinaryBladder",
  			       logicUrinaryBladder,
			       mother,
			       false,
			       0);

  // Sensitive Body Part
  if (sensitivity==true)
  { 
    G4SDManager* SDman = G4SDManager::GetSDMpointer();
    logicUrinaryBladder->SetSensitiveDetector( SDman->FindSensitiveDetector("BodyPartSD") );
  }

  // Visualization Attributes
  G4VisAttributes* UrinaryBladderVisAtt = new G4VisAttributes(G4Colour(0.85,0.65,0.125));
  UrinaryBladderVisAtt->SetForceSolid(true);
  logicUrinaryBladder->SetVisAttributes(UrinaryBladderVisAtt);

  G4cout << "UrinaryBladder created !!!!!!" << G4endl;

  // Testing UrinaryBladder Volume
  G4double UrinaryBladderVol = logicUrinaryBladder->GetSolid()->GetCubicVolume();
  G4cout << "Volume of UrinaryBladder = " << UrinaryBladderVol/cm3 << " cm^3" << G4endl;
  
  // Testing UrinaryBladder Material
  G4String UrinaryBladderMat = logicUrinaryBladder->GetMaterial()->GetName();
  G4cout << "Material of UrinaryBladder = " << UrinaryBladderMat << G4endl;
  
  // Testing Density
  G4double UrinaryBladderDensity = logicUrinaryBladder->GetMaterial()->GetDensity();
  G4cout << "Density of Material = " << UrinaryBladderDensity*cm3/g << " g/cm^3" << G4endl;

  // Testing Mass
  G4double UrinaryBladderMass = (UrinaryBladderVol)*UrinaryBladderDensity;
  G4cout << "Mass of UrinaryBladder = " << UrinaryBladderMass/gram << " g" << G4endl;

  
  return physUrinaryBladder;
}
