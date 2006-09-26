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
#include "G4MIRDBrain.hh"

#include "globals.hh"
#include "G4SDManager.hh"
#include "G4VisAttributes.hh"
#include "G4Ellipsoid.hh"
#include "G4ThreeVector.hh"
#include "G4VPhysicalVolume.hh"
#include "G4RotationMatrix.hh"
#include "G4Material.hh"
#include "G4LogicalVolume.hh"
#include "G4HumanPhantomMaterial.hh"
#include "G4VPhysicalVolume.hh"
#include "G4PVPlacement.hh"

G4MIRDBrain::G4MIRDBrain()
{
}

G4MIRDBrain::~G4MIRDBrain()
{

}

G4VPhysicalVolume* G4MIRDBrain::ConstructBrain(G4VPhysicalVolume* mother, G4String sex, G4bool sensitivity)
{

 G4HumanPhantomMaterial* material = new G4HumanPhantomMaterial();
 G4Material* soft = material -> GetMaterial("soft_tissue");
 delete material;

 G4double ax = 6.58 * cm;
 G4double by= 8.57 * cm;
 G4double cz = 5.73 * cm;

 G4Ellipsoid* brain = new G4Ellipsoid("Brain", ax, by, cz);
 

  G4LogicalVolume* logicBrain =  new G4LogicalVolume(brain, soft, 
						     "BrainVolume",
						     0, 0, 0);
  
  // Define rotation and position here!
  G4VPhysicalVolume* physBrain = new G4PVPlacement(0,
						   G4ThreeVector(0.*cm, 0.*cm, 8.25 * cm),
						   "physicalBrain",
						   logicBrain,
						   mother,
						   false,
						   0);
  // Sensitive Body Part
  if (sensitivity==true)
  { 
    G4SDManager* SDman = G4SDManager::GetSDMpointer();
    logicBrain->SetSensitiveDetector( SDman->FindSensitiveDetector("BodyPartSD") );
  }


  // Visualization Attributes
  G4VisAttributes* BrainVisAtt = new G4VisAttributes(G4Colour(0.41,0.41,0.41));
  BrainVisAtt->SetForceSolid(true);
  logicBrain->SetVisAttributes(BrainVisAtt);

  G4cout << "Brain created !!!!!!" << G4endl;

  // Testing Brain Volume
  G4double BrainVol = logicBrain->GetSolid()->GetCubicVolume();
  G4cout << "Volume of Brain = " << BrainVol/cm3 << " cm^3" << G4endl;
  
  // Testing Brain Material
  G4String BrainMat = logicBrain->GetMaterial()->GetName();
  G4cout << "Material of Brain = " << BrainMat << G4endl;
  
  // Testing Density
  G4double BrainDensity = logicBrain->GetMaterial()->GetDensity();
  G4cout << "Density of Material = " << BrainDensity*cm3/g << " g/cm^3" << G4endl;

  // Testing Mass
  G4double BrainMass = (BrainVol)*BrainDensity;
  G4cout << "Mass of Brain = " << BrainMass/gram << " g" << G4endl;

  
  return physBrain;
}
