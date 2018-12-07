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
//
#include "G4MIRDPelvis.hh"

#include "globals.hh"
#include "G4SystemOfUnits.hh"
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
#include "G4HumanPhantomColour.hh"

G4MIRDPelvis::G4MIRDPelvis()
{
}

G4MIRDPelvis::~G4MIRDPelvis()
{

}


G4VPhysicalVolume* G4MIRDPelvis::Construct(const G4String& volumeName,G4VPhysicalVolume* mother, 
					   const G4String& colourName, G4bool wireFrame,G4bool)
{
  G4HumanPhantomMaterial* material = new G4HumanPhantomMaterial();
   
  G4cout<<"Construct "<<volumeName<<" with mother volume "<<mother->GetName()<<G4endl;
   
  G4Material* skeleton = material -> GetMaterial("skeleton");
 
  delete material;
  /*
    G4double dx= 10.35 * cm;//12. *cm; // a2 
    G4double dy= 11.76 * cm;//12. * cm; // b2
    G4double dz= 9.915 * cm; // z2/2

    G4VSolid* outPelvis = new G4EllipticalTube("OutPelvis",dx, dy, dz);

    G4double dx_in = 9.75 * cm;//11.3 * cm; // a1
    G4double dy_in = 11.07 *cm; //11.3* cm; //b1
    G4double dz_in = 10. * cm;//11.0 *cm; // z2/2

  */
  G4double dx= 12. *cm; // a2
  G4double dy= 12. * cm; //b2
  G4double dz= 11. * cm; // z2/2

  G4VSolid* outPelvis = new G4EllipticalTube("OutPelvis",dx, dy, dz);

  dx = 11.3 * cm; // a1
  dy = 11.3* cm; // b1
  dz = 12.0 *cm; // z2/2
 
  G4VSolid* inPelvis = new G4EllipticalTube("InPelvis",dx, dy, dz);

  G4double x = 28. * cm; // a2 * 2
  G4double y = 28. * cm; //b2*2
  G4double z = 24. *cm; // z2

  G4VSolid* subPelvis = new G4Box("SubtrPelvis", x/2., y/2., z/2.);

  

  G4SubtractionSolid* firstPelvis = new G4SubtractionSolid("FirstPelvis",
							   outPelvis,
							   inPelvis, 0, G4ThreeVector(0.*cm, -0.8 *cm, 0. * cm)); 
							   
    
  G4SubtractionSolid* secondPelvis = new G4SubtractionSolid("SecondPelvis",
							    firstPelvis,
							    subPelvis, 0, 
							    G4ThreeVector(0.0,
									  -14. * cm, 0.*cm));
  // half of the y size of the box
 
  
  G4SubtractionSolid* pelvis = new G4SubtractionSolid("Pelvis", secondPelvis, subPelvis,
						      0, 
						      G4ThreeVector(0.0,
								    22. * cm, -9. *cm)); 

 
  G4LogicalVolume* logicPelvis = new G4LogicalVolume(pelvis, skeleton,
						     "logical" + volumeName, 0, 0, 0);
  
 
  G4VPhysicalVolume* physPelvis = new G4PVPlacement(0,G4ThreeVector(0.0, -3. * cm,-24. * cm),// 0, y02, z position
						    // with respect to the trunk 
						    "physicalPelvis",
						    logicPelvis,
						    mother,
						    false,
						    0, true);

 

  // Visualization Attributes
  //  G4VisAttributes* PelvisVisAtt = new G4VisAttributes(G4Colour(0.46,0.53,0.6));
 
  G4HumanPhantomColour* colourPointer = new G4HumanPhantomColour();
  G4Colour colour = colourPointer -> GetColour(colourName);
  G4VisAttributes* PelvisVisAtt = new G4VisAttributes(colour);
  PelvisVisAtt->SetForceSolid(wireFrame);
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
