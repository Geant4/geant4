//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
//
// Previous authors: G. Guerrieri, S. Guatelli and M. G. Pia, INFN Genova, Italy
// Authors (since 2007): S. Guatelli, University of Wollongong, Australia
// 

#include "G4MIRDRightLung.hh"
#include "globals.hh"
#include "G4SystemOfUnits.hh"
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
#include "G4Box.hh"
#include "G4HumanPhantomColour.hh"

G4VPhysicalVolume* G4MIRDRightLung::Construct(const G4String& volumeName,G4VPhysicalVolume* mother, 
					      const G4String& colourName, G4bool wireFrame,G4bool)
{
  G4cout << "Construct " << volumeName << " with mother volume " <<mother->GetName()<<G4endl;
 
  auto* material = new G4HumanPhantomMaterial();
  auto* lung_material = material -> GetMaterial("lung_material");
  delete material;

  G4double ax = 5. *cm; //a
  G4double by = 7.5 *cm; //b
  G4double cz = 24.*cm; //c
  G4double zcut1 = 0.0 *cm; 
  G4double zcut2=24. *cm;
 
  auto* oneLung = new G4Ellipsoid("OneLung",ax, by, cz, zcut1,zcut2);
 
  ax= 5.*cm; 
  by= 7.5*cm; 
  cz= 24.*cm;


  auto* subtrLung = new G4Ellipsoid("subtrLung",ax, by, cz);

  // y<0

  G4double dx = 5.5* cm;
  G4double dy = 8.5 * cm;
  G4double dz = 24. * cm;

  auto* box = new G4Box("Box", dx, dy, dz);
 
  auto* section = new G4SubtractionSolid("BoxSub", subtrLung, box, nullptr, G4ThreeVector(0.*cm, 8.5* cm, 0.*cm)); 
  //G4SubtractionSolid* section2 = new G4SubtractionSolid("BoxSub2", subtrLung, box, 0, G4ThreeVector(0.*cm, -8.5* cm, 0.*cm)); 
 
  auto* lung1 =  new G4SubtractionSolid("Lung1", oneLung,
					  section,
				          nullptr, G4ThreeVector(6.*cm,0*cm,0.0*cm));
 
  auto* logicRightLung = new G4LogicalVolume(lung1,lung_material,
							"logical" + volumeName, nullptr, nullptr, nullptr); 
 
  auto* rm = new G4RotationMatrix();
  rm->rotateZ(-360.*degree);
  G4VPhysicalVolume* physRightLung = new G4PVPlacement(rm,G4ThreeVector(-8.50 *cm, 0.0*cm, 8.5*cm),
						       "physicalRightLung",                    
						       logicRightLung,
						       mother,
						       false,
						       0, true);

  // Visualization Attributes
  // G4VisAttributes* RightLungVisAtt = new G4VisAttributes(G4Colour(0.25,0.41,0.88));
  auto* colourPointer = new G4HumanPhantomColour();
  G4Colour colour = colourPointer -> GetColour(colourName);
  auto* RightLungVisAtt = new G4VisAttributes(colour);
  RightLungVisAtt->SetForceSolid(wireFrame);
  logicRightLung->SetVisAttributes(RightLungVisAtt); 

  G4cout << "RightLung created !!!!!!" << G4endl;

  // Testing RightLung Volume
  G4double RightLungVol = logicRightLung->GetSolid()->GetCubicVolume();
 
  G4cout << "Volume of RightLung = " << (RightLungVol)/cm3 << " cm^3" << G4endl;
  
  // Testing RightLung Material
  G4String RightLungMat = logicRightLung->GetMaterial()->GetName();
  G4cout << "Material of RightLung = " << RightLungMat << G4endl;
  
  // Testing Density
  G4double RightLungDensity = logicRightLung->GetMaterial()->GetDensity();
  G4cout << "Density of Material = " << RightLungDensity*cm3/g << " g/cm^3" << G4endl;

  // Testing Mass
  G4double RightLungMass = (RightLungVol)*RightLungDensity;
  G4cout << "Mass of RightLung = " << RightLungMass/gram << " g" << G4endl;
  
  return physRightLung;
}
