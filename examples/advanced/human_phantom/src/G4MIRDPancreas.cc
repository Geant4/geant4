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
#include "G4MIRDPancreas.hh"

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
#include "G4Box.hh"
#include "G4SubtractionSolid.hh"

G4MIRDPancreas::G4MIRDPancreas()
{
}

G4MIRDPancreas::~G4MIRDPancreas()
{

}

G4VPhysicalVolume* G4MIRDPancreas::ConstructPancreas(G4VPhysicalVolume* mother, G4String sex, G4bool sensitivity)
{

  G4cout << "ConstructPancreas for " << sex << G4endl;
 
 G4HumanPhantomMaterial* material = new G4HumanPhantomMaterial();
 G4Material* soft = material -> GetMaterial("soft_tissue");
 delete material;

  G4double ax= 2.87*cm;
  G4double by= 1.14*cm;
  G4double cz= 13.32*cm;
  G4double zcut1= -13.32 *cm;
  G4double zcut2= 0.0 *cm; 

  G4Ellipsoid* pancreasFirst =  new G4Ellipsoid("PancreasFirst",ax, by, cz,
						zcut1, zcut2);

  G4double xx = 5.74 * cm;
  G4double yy = 2.30 * cm;
  G4double zz = 10.71 * cm;
  G4Box* subtrPancreas = new G4Box("SubtrPancreas",xx/2., yy/2., zz/2.);

  G4SubtractionSolid* pancreas = new G4SubtractionSolid("pancreas",
							pancreasFirst,
							subtrPancreas,
							0, 
							G4ThreeVector(-2.87 * cm,0.0,-8.685*cm));
 
  G4LogicalVolume* logicPancreas = new G4LogicalVolume(pancreas, soft,
						       "PancreasVolume",
						       0, 0, 0);
  G4RotationMatrix* rotation = new G4RotationMatrix();
  rotation ->rotateY(90. * degree);


  G4VPhysicalVolume* physPancreas = new G4PVPlacement(rotation,
						      G4ThreeVector(-0.72 *cm, 0.0, 1.8*cm),
      			       "physicalPancreas",
  			       logicPancreas,
			       mother,
			       false,
			       0);
  // Sensitive Body Part
  if (sensitivity==true)
  { 
    G4SDManager* SDman = G4SDManager::GetSDMpointer();
    logicPancreas->SetSensitiveDetector( SDman->FindSensitiveDetector("BodyPartSD") );
  }

  // Visualization Attributes
  G4VisAttributes* PancreasVisAtt = new G4VisAttributes(G4Colour(0.28,0.82,0.8));
  PancreasVisAtt->SetForceSolid(true);
  logicPancreas->SetVisAttributes(PancreasVisAtt);

  G4cout << "Pancreas created !!!!!!" << G4endl;

  // Testing Pancreas Volume
  G4double PancreasVol = logicPancreas->GetSolid()->GetCubicVolume();
  G4cout << "Volume of Pancreas = " << PancreasVol/cm3 << " cm^3" << G4endl;
  
  // Testing Pancreas Material
  G4String PancreasMat = logicPancreas->GetMaterial()->GetName();
  G4cout << "Material of Pancreas = " << PancreasMat << G4endl;
  
  // Testing Density
  G4double PancreasDensity = logicPancreas->GetMaterial()->GetDensity();
  G4cout << "Density of Material = " << PancreasDensity*cm3/g << " g/cm^3" << G4endl;

  // Testing Mass
  G4double PancreasMass = (PancreasVol)*PancreasDensity;
  G4cout << "Mass of Pancreas = " << PancreasMass/gram << " g" << G4endl;

  
  return physPancreas;
}
