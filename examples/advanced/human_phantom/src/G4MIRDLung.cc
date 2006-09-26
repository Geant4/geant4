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
#include "G4UnionSolid.hh"
#include "G4SubtractionSolid.hh"

G4MIRDLung::G4MIRDLung()
{
}

G4MIRDLung::~G4MIRDLung()
{

}

G4VPhysicalVolume* G4MIRDLung::ConstructLung(G4VPhysicalVolume* mother, G4String sex, G4bool sensitivity)
{

 G4cout << "ConstructLungs for " << sex << G4endl;
 
 G4HumanPhantomMaterial* material = new G4HumanPhantomMaterial();
 G4Material* lung_material = material -> GetMaterial("lung_material");
 delete material;

 G4double ax = 4.09 *cm;
 G4double by = 6.98 *cm; 
 G4double cz = 20.55*cm;
 G4double zcut1 = 0.0 *cm; 
 G4double zcut2=20.55 *cm;
 
 G4Ellipsoid* oneLung = new G4Ellipsoid("OneLung",ax, by, cz, zcut1,zcut2);
 
 ax= 5.0*cm;
 by= 7.0*cm;
 cz= 20.55*cm;

 G4Ellipsoid* subtrLung = new G4Ellipsoid("subtrLung",ax, by, cz);

 G4UnionSolid* lung1 = new G4UnionSolid("Lung1", oneLung, oneLung,0,
				    G4ThreeVector(14.66*cm, 0.0*cm,0.*cm));

 G4SubtractionSolid* lung =  new G4SubtractionSolid("Lung", lung1,
					       subtrLung,
					       0, G4ThreeVector(7.50*cm,-4.0*cm,0.0));

 G4LogicalVolume* logicLung = new G4LogicalVolume(lung,lung_material,
						  "LungVolume", 0, 0, 0); 
  

  G4VPhysicalVolume* physLung = new G4PVPlacement(0,G4ThreeVector(-7.33 *cm, 0.0*cm, 7.66*cm),
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
