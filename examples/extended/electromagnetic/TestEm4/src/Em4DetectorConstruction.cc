// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: Em4DetectorConstruction.cc,v 1.2 1999-12-15 14:49:06 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#include "Em4DetectorConstruction.hh"

#include "G4Material.hh"
#include "G4Tubs.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"

#include "G4VisAttributes.hh"
#include "G4Colour.hh"

#include "G4ios.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

Em4DetectorConstruction::Em4DetectorConstruction()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

Em4DetectorConstruction::~Em4DetectorConstruction()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4VPhysicalVolume* Em4DetectorConstruction::Construct()
{
  //
  // define a material from elements.   case 1: chemical molecule
  //
 
  G4String name, symbol;             //a=mass of a mole;
  G4double a, z;                     //z=mean number of protons;  
  G4int ncomponents, natoms;
 
  a = 12.01*g/mole;
  G4Element* elC  = new G4Element(name="Carbon"  ,symbol="C" , z= 6., a);

  a = 18.99*g/mole;
  G4Element* elF  = new G4Element(name="Fluorine",symbol="N" , z= 9., a);
 
  G4double density = 1.61*g/cm3;
  G4Material* matC6F6 = new G4Material(name="FluorCarbonate",density,ncomponents=2);
  matC6F6->AddElement(elC, natoms=6);
  matC6F6->AddElement(elF, natoms=6);
  
  G4cout << matC6F6 << G4endl;
  
  //     
  // Container
  //
  
  G4double Rmin=0., Rmax=5*cm, deltaZ= 5*cm, Phimin=0., deltaPhi=360*degree;

  G4Tubs*  
  solidWorld = new G4Tubs("C6F6",			//its name
                   Rmin,Rmax,deltaZ,Phimin,deltaPhi);	//its size

  G4LogicalVolume*                         
  logicWorld = new G4LogicalVolume(solidWorld,		//its solid
                                   matC6F6,		//its material
                                   "C6F6");		//its name
  G4VPhysicalVolume*                                   
  physiWorld = new G4PVPlacement(0,			//no rotation
  				 G4ThreeVector(),	//at (0,0,0)
                                 "C6F6",		//its name
                                 logicWorld,		//its logical volume
                                 NULL,			//its mother  volume
                                 false,			//no boolean operation
                                 0);			//copy number
     
  //                                        
  // Visualization attributes
  //
  
  G4VisAttributes* visAtt= new G4VisAttributes(G4Colour(1.0,1.0,1.0));
  visAtt->SetVisibility(true);
  logicWorld->SetVisAttributes(visAtt);
  
  //
  //always return the physical World
  //
  
  return physiWorld;
}
