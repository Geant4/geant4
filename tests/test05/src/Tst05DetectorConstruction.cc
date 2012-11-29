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
// $Id$
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo...... 
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
 
#include "Tst05DetectorConstruction.hh"

#include "G4Material.hh"
#include "G4Box.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4PVParameterised.hh"
#include "G4SDManager.hh"

#include "G4UserLimits.hh"

#include "G4VisAttributes.hh"
#include "G4Colour.hh"

#include "G4ios.hh"
#include "G4SystemOfUnits.hh"

#include "G4TriangularFacet.hh"
#include "G4QuadrangularFacet.hh"


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
 
Tst05DetectorConstruction::Tst05DetectorConstruction()
:solidWorld(0),  logicWorld(0),  physiWorld(0),
 solidTarget(0), logicTarget(0), physiTarget(0)
{;}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
 
Tst05DetectorConstruction::~Tst05DetectorConstruction()
{;}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
 
G4VPhysicalVolume* Tst05DetectorConstruction::Construct()
{
//--------- Material definition ---------

  G4double a, z;
  G4double density;
  G4int nel;

  //Air
  G4Element* N = new G4Element("Nitrogen", "N", z=7., a= 14.01*g/mole);
  G4Element* O = new G4Element("Oxygen"  , "O", z=8., a= 16.00*g/mole);
   
  G4Material* Air = new G4Material("Air", density= 1.29*mg/cm3, nel=2);
  Air->AddElement(N, 70*perCent);
  Air->AddElement(O, 30*perCent);

  //Lead
  G4Material* Pb = 
  new G4Material("Lead", z=82., a= 207.19*g/mole, density= 11.35*g/cm3);
    
  // Print all the materials defined.
  //
  G4cout << G4endl << "The materials defined are : " << G4endl << G4endl;
  G4cout << *(G4Material::GetMaterialTable()) << G4endl;

//--------- Sizes of the principal geometrical components (solids)  ---------
  
  fTargetLength  = 10.0 * cm;                        // Full length of Target
  
  fWorldLength= 1.2 *(fTargetLength*2.0);
   
  G4double targetSize  = 0.5*fTargetLength;    // Half length of the Target  
      
//--------- Definitions of Solids, Logical Volumes, Physical Volumes ---------
  
  //------------------------------ 
  // World
  //------------------------------ 

  G4double HalfWorldLength = 0.5*fWorldLength;
  
 solidWorld= new G4Box("world",HalfWorldLength,HalfWorldLength,HalfWorldLength);
 logicWorld= new G4LogicalVolume( solidWorld, Air, "World", 0, 0, 0);
  
  //  Must place the World Physical volume unrotated at (0,0,0).
  // 
  physiWorld = new G4PVPlacement(0,               // no rotation
                                 G4ThreeVector(), // at (0,0,0)
                                 logicWorld,      // its logical volume
				 "World",         // its name
                                 0,               // its mother  volume
                                 false,           // no boolean operations
                                 0);              // no field specific to volume
				 
  //------------------------------ 
  // Target
  //------------------------------
  
  G4ThreeVector positionTarget = G4ThreeVector(0.0,0.0,0.0);
   
  //solidTarget = new G4Box("target",targetSize,targetSize,targetSize);

  solidTarget = new G4TessellatedSolid("target");
  G4TriangularFacet *facet1 = new
  G4TriangularFacet (G4ThreeVector(-targetSize,-targetSize,        0.0),
                     G4ThreeVector(+targetSize,-targetSize,        0.0),
                     G4ThreeVector(        0.0,        0.0,+targetSize),
                     ABSOLUTE);
  G4TriangularFacet *facet2 = new
  G4TriangularFacet (G4ThreeVector(+targetSize,-targetSize,        0.0),
                     G4ThreeVector(+targetSize,+targetSize,        0.0),
                     G4ThreeVector(        0.0,        0.0,+targetSize),
                     ABSOLUTE);
  G4TriangularFacet *facet3 = new
  G4TriangularFacet (G4ThreeVector(+targetSize,+targetSize,        0.0),
                     G4ThreeVector(-targetSize,+targetSize,        0.0),
                     G4ThreeVector(        0.0,        0.0,+targetSize),
                     ABSOLUTE);
  G4TriangularFacet *facet4 = new
  G4TriangularFacet (G4ThreeVector(-targetSize,+targetSize,        0.0),
                     G4ThreeVector(-targetSize,-targetSize,        0.0),
                     G4ThreeVector(        0.0,        0.0,+targetSize),
                     ABSOLUTE);
  G4QuadrangularFacet *facet5 = new
  G4QuadrangularFacet (G4ThreeVector(-targetSize,-targetSize,        0.0),
                     G4ThreeVector(-targetSize,+targetSize,        0.0),
                     G4ThreeVector(+targetSize,+targetSize,        0.0),
                     G4ThreeVector(+targetSize,-targetSize,        0.0),
                     ABSOLUTE);

  solidTarget->AddFacet((G4VFacet*) facet1);
  solidTarget->AddFacet((G4VFacet*) facet2);
  solidTarget->AddFacet((G4VFacet*) facet3);
  solidTarget->AddFacet((G4VFacet*) facet4);
  solidTarget->AddFacet((G4VFacet*) facet5);
  
  solidTarget->SetSolidClosed(true);
  logicTarget = new G4LogicalVolume(solidTarget,Pb,"Target",0,0,0);
  physiTarget = new G4PVPlacement(0,               // no rotation
				  positionTarget,  // at (x,y,z)
				  logicTarget,     // its logical volume				  
				  "Target",        // its name
				  logicWorld,      // its mother  volume
				  false,           // no boolean operations
				  0);              // no particular field 
  
//--------- Visualization attributes -------------------------------

  G4VisAttributes* BoxVisAtt= new G4VisAttributes(G4Colour(1.0,1.0,1.0));
  BoxVisAtt->SetVisibility(false);
  G4VisAttributes* FacetVisAtt= new G4VisAttributes(G4Colour(1.0,0.0,0.0));
  FacetVisAtt->SetVisibility(true);
  
  logicWorld  ->SetVisAttributes(BoxVisAtt);  
  logicTarget ->SetVisAttributes(FacetVisAtt);
  
//--------- example of User Limits -------------------------------

  // below is an example of how to set tracking constraints in a given
  // logical volume(see also in N02PhysicsList how to setup the process
  // G4UserSpecialCuts).  
  // Sets a max Step length in the tracker region
  // G4double maxStep = 0.5*ChamberWidth, maxLength = 2*fTrackerLength;
  // G4double maxTime = 0.1*ns, minEkin = 10*MeV;
  // logicTracker->SetUserLimits(new G4UserLimits(maxStep,maxLength,maxTime,
  //                                               minEkin));
  
  return physiWorld;
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
