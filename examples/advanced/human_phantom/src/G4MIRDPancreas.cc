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
#include "G4MIRDPancreas.hh"

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
#include "G4Box.hh"
#include "G4SubtractionSolid.hh"
#include "G4HumanPhantomColour.hh"

G4VPhysicalVolume* G4MIRDPancreas::Construct(const G4String& volumeName,G4VPhysicalVolume* mother,
					     const G4String& colourName
					     ,G4bool wireFrame, G4bool)
{
  G4cout<<"Construct "<<volumeName<<" with mother volume "<<mother->GetName()<<G4endl;
  
  auto* material = new G4HumanPhantomMaterial();
  auto* soft = material -> GetMaterial("soft_tissue");
  delete material;
  
  G4double ax= 3.*cm; //c
  G4double by= 1.*cm;//b
  G4double cz= 15.*cm;//a
  G4double zcut1= -15. *cm;// -a
  G4double zcut2= 0.0 *cm; 
  
  auto* pancreasFirst =  new G4Ellipsoid("PancreasFirst",ax, by, cz,
						zcut1, zcut2);
  
  G4double xx = 6. * cm;// 2*c
  G4double yy = 2. * cm;// 2*b
  G4double zz = 12. * cm; // cz - x1 = 3 cm
  auto* subtrPancreas = new G4Box("SubtrPancreas",xx/2., yy/2., zz/2.);
  
  auto* pancreas = new G4SubtractionSolid("pancreas",
							pancreasFirst,
							subtrPancreas,
							nullptr, 
							G4ThreeVector(-3 * cm,0.0,-9.*cm));
  
  auto* logicPancreas = new G4LogicalVolume(pancreas, soft,
						       "logical" + volumeName,
						         nullptr, nullptr, nullptr);
  
  auto* rm = new G4RotationMatrix();
  rm->rotateY(90.*degree);
  G4VPhysicalVolume* physPancreas = new G4PVPlacement(rm,
						      G4ThreeVector(-0. *cm, 0.0, 2*cm),//x0, 0, 2 cm
						      "physicalPancreas",
						      logicPancreas,
						      mother,
						      false,
						      0, true);
 
  
  // Visualization Attributes
  // G4VisAttributes* PancreasVisAtt = new G4VisAttributes(G4Colour(0.28,0.82,0.8));
  auto* colourPointer = new G4HumanPhantomColour();
  G4Colour colour = colourPointer -> GetColour(colourName);
  auto* PancreasVisAtt = new G4VisAttributes(colour);
  PancreasVisAtt->SetForceSolid(wireFrame);
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
