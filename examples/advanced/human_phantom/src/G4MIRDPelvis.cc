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
// Authors: S. Guatelli and M. G. Pia, INFN Genova, Italy
// 
// Based on code developed by the undergraduate student G. Guerrieri 
// Note: this is a preliminary beta-version of the code; an improved 
// version will be distributed in the next Geant4 public release, compliant
// with the design in a forthcoming publication, and subject to a 
// design and code review.
//
#include "G4MIRDPelvis.hh"

#include "globals.hh"
#include "G4SDManager.hh"
#include "G4VisAttributes.hh"
#include "G4HumanPhantomMaterial.hh"
#include "G4EllipticalTube.hh"
#include "G4RotationMatrix.hh"
#include "G4ThreeVector.hh"
#include "G4VPhysicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4SubtractionSolid.hh"
#include "G4Box.hh"
#include "G4VSolid.hh"
#include "G4LogicalVolume.hh"

G4MIRDPelvis::G4MIRDPelvis()
{
}

G4MIRDPelvis::~G4MIRDPelvis()
{

}

G4VPhysicalVolume* G4MIRDPelvis::ConstructPelvis(G4VPhysicalVolume* mother, G4String sex, G4bool sensitivity)
{
   G4HumanPhantomMaterial* material = new G4HumanPhantomMaterial();
   
  G4cout << "ConstructMiddleLowerSpine for "<< sex <<G4endl;
   
  G4Material* skeleton = material -> GetMaterial("skeleton");
 
  delete material;

  G4double dx= 10.35 *cm;
  G4double dy= 11.76 * cm;
  G4double dz= 9.915 * cm;

  G4VSolid* outPelvis = new G4EllipticalTube("OutPelvis",dx, dy, dz);

  dx = 9.75 * cm;
  dy = 11.07* cm;
  dz = 10.0 *cm;
 
  G4VSolid* inPelvis = new G4EllipticalTube("InPelvis",dx, dy, dz);

  G4double x = 20.71 * cm;
  G4double y = 23.52 * cm;
  G4double z = 19.84 *cm;

  G4VSolid* subPelvis = new G4Box("SubtrPelvis", x/2., y/2., z/2.);

  G4SubtractionSolid* firstPelvis = new G4SubtractionSolid("FirstPelvis",
							   outPelvis,
							   inPelvis); 
							   
 
  G4SubtractionSolid* secondPelvis = new G4SubtractionSolid("SecondPelvis",
							    firstPelvis,
							    subPelvis, 0, 
							    G4ThreeVector(0.0,
									  14.70 * cm, -8.21 *cm));

  G4SubtractionSolid* pelvis = new G4SubtractionSolid("Pelvis", secondPelvis, subPelvis,
						      0, 
						      G4ThreeVector(0.0,
								    -11.76 * cm, 0.0));


  G4LogicalVolume* logicPelvis = new G4LogicalVolume(pelvis, skeleton,
						       "PelvisVolume", 0, 0, 0);
  
 
  G4VPhysicalVolume* physPelvis = new G4PVPlacement(0,G4ThreeVector(0.0, -2.94 * cm,-21.635 * cm),
      			       "physicalPelvis",
  			       logicPelvis,
			       mother,
			       false,
			       0);

  // Sensitive Body Part
  if (sensitivity==true)
  { 
    G4SDManager* SDman = G4SDManager::GetSDMpointer();
    logicPelvis->SetSensitiveDetector( SDman->FindSensitiveDetector("BodyPartSD") );
  }

  // Visualization Attributes
  G4VisAttributes* PelvisVisAtt = new G4VisAttributes(G4Colour(0.46,0.53,0.6));
  PelvisVisAtt->SetForceSolid(false);
  logicPelvis->SetVisAttributes(PelvisVisAtt);

  G4cout << "Pelvis created !!!!!!" << G4endl;

  // Testing Pelvis Volume
  G4double PelvisVol = logicPelvis->GetSolid()->GetCubicVolume();
  G4cout << "Volume of Pelvis = " << PelvisVol/cm3 << " cm^3" << G4endl;
  
  // Testing Pelvis Material
  G4String PelvisMat = logicPelvis->GetMaterial()->GetName();
  G4cout << "Material of Pelvis = " << PelvisMat << G4endl;
  
  // Testing Density
  G4double PelvisDensity = logicPelvis->GetMaterial()->GetDensity();
  G4cout << "Density of Material = " << PelvisDensity*cm3/g << " g/cm^3" << G4endl;

  // Testing Mass
  G4double PelvisMass = (PelvisVol)*PelvisDensity;
  G4cout << "Mass of Pelvis = " << PelvisMass/gram << " g" << G4endl;

  
  return physPelvis;
}
