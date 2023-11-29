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
// Previous authors: G. Guerrieri, S. Guatelli and M. G. Pia, INFN Genova, Italy
// Authors (since 2007): S. Guatelli, University of Wollongong, Australia
// 
//
#include "G4MIRDLeftLung.hh"

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

G4VPhysicalVolume* G4MIRDLeftLung::Construct(const G4String& volumeName,G4VPhysicalVolume* mother, 
					     const G4String& colourName, G4bool wireFrame,G4bool)
{

  G4cout<<"Construct "<<volumeName<<" with mother volume "<<mother->GetName()<<G4endl;
 
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
 

  //G4SubtractionSolid* section = new G4SubtractionSolid("BoxSub", subtrLung, box, 0, G4ThreeVector(0.*cm, 8.5* cm, 0.*cm)); 
  auto* section2 = new G4SubtractionSolid("BoxSub2", subtrLung, box, nullptr, G4ThreeVector(0.*cm, 8.5* cm, 0.*cm)); 

  //G4SubtractionSolid* lung1 =  new G4SubtractionSolid("Lung1", oneLung,
  //				       section,
  //				       0, G4ThreeVector(6.*cm,0*cm,0.0*cm));
 
  auto* lung2 =  new G4SubtractionSolid("Lung2", oneLung,
						      section2,
						      nullptr, G4ThreeVector(-6.*cm,0*cm,0.0*cm));

  // G4RotationMatrix* matrix = new G4RotationMatrix();  
  // matrix->rotateX(180. * degree);
  //matrix ->rotateZ(180.*degree);
  //matrix -> rotateY(180.* degree);

  //G4UnionSolid* lungs = new G4UnionSolid("Lungs", lung1, lung2, matrix, G4ThreeVector(17*cm, 0., 0.));


  auto* logicLeftLung = new G4LogicalVolume(lung2,lung_material,
						       "logical" + volumeName, nullptr, nullptr, nullptr); 
  

  G4VPhysicalVolume* physLeftLung = new G4PVPlacement(nullptr,G4ThreeVector(8.50 *cm, 0.0*cm, 8.5*cm),
						      "physicalLeftLung",                    
						      logicLeftLung,
						      mother,
						      false,
						      0, true);


  // Visualization Attributes
  //G4VisAttributes* LeftLungVisAtt = new G4VisAttributes(G4Colour(0.25,0.41,0.88));
  auto* colourPointer = new G4HumanPhantomColour();
  G4Colour colour = colourPointer -> GetColour(colourName);
  auto* LeftLungVisAtt = new G4VisAttributes(colour);
  LeftLungVisAtt->SetForceSolid(wireFrame);
  logicLeftLung->SetVisAttributes(LeftLungVisAtt); 

  G4cout << "LeftLung created !!!!!!" << G4endl;

  // Testing LeftLung Volume
  G4double LeftLungVol = logicLeftLung->GetSolid()->GetCubicVolume();
 
  G4cout << "Volume of LeftLung = " << (LeftLungVol)/cm3 << " cm^3" << G4endl;
  
  // Testing LeftLung Material
  G4String LeftLungMat = logicLeftLung->GetMaterial()->GetName();
  G4cout << "Material of LeftLung = " << LeftLungMat << G4endl;
  
  // Testing Density
  G4double LeftLungDensity = logicLeftLung->GetMaterial()->GetDensity();
  G4cout << "Density of Material = " << LeftLungDensity*cm3/g << " g/cm^3" << G4endl;

  // Testing Mass
  G4double LeftLungMass = (LeftLungVol)*LeftLungDensity;
  G4cout << "Mass of LeftLung = " << LeftLungMass/gram << " g" << G4endl;
  
  return physLeftLung;
}
