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
#include "G4MIRDMaleGenitalia.hh"

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
#include "G4HumanPhantomColour.hh"
#include "G4Trap.hh"
#include "G4Cons.hh"
#include "G4SubtractionSolid.hh"

G4MIRDMaleGenitalia::G4MIRDMaleGenitalia()
{
}

G4MIRDMaleGenitalia::~G4MIRDMaleGenitalia()
{

}


G4VPhysicalVolume* G4MIRDMaleGenitalia::Construct(const G4String& volumeName,G4VPhysicalVolume* mother, 
						  const G4String& colourName, G4bool wireFrame, G4bool)
{ 
  G4cout<<"Construct "<<volumeName<<" with mother volume "<<mother->GetName()<<G4endl;

  G4HumanPhantomMaterial* material = new G4HumanPhantomMaterial();
  G4Material* soft = material -> GetMaterial("soft_tissue");
  delete material;
 

  G4double pDz=2.4*cm;
  G4double pTheta=0*degree;
  G4double pPhi=0*degree;
  G4double pDy1=4.76*cm;
  G4double pDx1=9.52*cm; 
  G4double pDx2=9.52*cm;
  G4double pAlp1=0*degree;
  G4double pDy2=5*cm;
  G4double pDx3=10*cm; 
  G4double pDx4=10*cm;
  G4double pAlp2=0*degree;

  G4Trap* genitaliaTrap= new G4Trap("GenitaliaTrap", 
				    pDz,pTheta,pPhi,pDy1,
				    pDx1,pDx2,pAlp1,pDy2,
				    pDx3,pDx4,pAlp2);
 

  G4double rmin1 = 0.* cm;
  G4double rmin2 = 0.* cm;
  G4double dz= 5 * cm; 
  G4double rmax1= 9.51 * cm;
  G4double rmax2= 10.01 * cm;
  G4double startphi= 0.* degree;
  G4double deltaphi= 360. * degree;

  G4Cons* genitaliaLegL = new G4Cons("GenitaliaLegL",  
				     rmin1, rmax1, 
				     rmin2, rmax2, dz/2., 
				     startphi, deltaphi);

  G4Cons* genitaliaLegR = new G4Cons("GenitaliaLegR",  
				     rmin1, rmax1, 
				     rmin2, rmax2, dz/2., 
				     startphi, deltaphi);
 
  G4UnionSolid* genitaliaLegs = new G4UnionSolid("GenitaliaLegs",genitaliaLegL,genitaliaLegR,
						 0,
						 G4ThreeVector(20.* cm, 0.*cm,0* cm) );



  G4SubtractionSolid* MaleGenitalia = new G4SubtractionSolid("MaleGenitalia",genitaliaTrap,genitaliaLegs,
							     0,//
							     G4ThreeVector(-10.* cm, -5.*cm,0* cm) );
 
 
  G4LogicalVolume* logicMaleGenitalia = new G4LogicalVolume(MaleGenitalia,
							    soft,
							    "logical" + volumeName,
							    0, 0, 0);
 
  // Define rotation and position here!
  G4VPhysicalVolume* physMaleGenitalia = new G4PVPlacement(0,
							   G4ThreeVector(0*cm,5.*cm, -2.4*cm),
							   "physicalMaleGenitalia",
							   logicMaleGenitalia,
							   mother,
							   false,
							   0, true);
 
  // Visualization Attributes
  //G4VisAttributes* MaleGenitaliaVisAtt = new G4VisAttributes(G4Colour(0.85,0.44,0.84));
  G4HumanPhantomColour* colourPointer = new G4HumanPhantomColour();
  G4Colour colour = colourPointer -> GetColour(colourName);
  delete colourPointer;
 
  G4VisAttributes* MaleGenitaliaVisAtt = new G4VisAttributes(colour);
  MaleGenitaliaVisAtt->SetForceSolid(wireFrame);
  logicMaleGenitalia->SetVisAttributes(MaleGenitaliaVisAtt);
 
  G4cout << "MaleGenitalia created !!!!!!" << G4endl;
 
  // Testing MaleGenitalia Volume
  G4double MaleGenitaliaVol = logicMaleGenitalia->GetSolid()->GetCubicVolume();
  G4cout << "Volume of MaleGenitalia = " << MaleGenitaliaVol/cm3 << " cm^3" << G4endl;
 
  // Testing MaleGenitalia Material
  G4String MaleGenitaliaMat = logicMaleGenitalia->GetMaterial()->GetName();
  G4cout << "Material of MaleGenitalia = " << MaleGenitaliaMat << G4endl;
 
  // Testing Density
  G4double MaleGenitaliaDensity = logicMaleGenitalia->GetMaterial()->GetDensity();
  G4cout << "Density of Material = " << MaleGenitaliaDensity*cm3/g << " g/cm^3" << G4endl;
 
  // Testing Mass
  G4double MaleGenitaliaMass = (MaleGenitaliaVol)*MaleGenitaliaDensity;
  G4cout << "Mass of MaleGenitalia = " << MaleGenitaliaMass/gram << " g" << G4endl;
 
  return physMaleGenitalia;
}
