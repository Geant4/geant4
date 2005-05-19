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
//
// $Id: SCDetectorConstruction.cc,v 1.1 2005-05-19 13:07:29 link Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo...... 
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
 
#include "SCDetectorConstruction.hh"
#include "SCDetectorMessenger.hh"
#include "SCMagneticField.hh"
#include "SCTrackerSD.hh"

#include "G4Material.hh"
#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4Torus.hh"
#include "G4TwistedTrap.hh"
#include "G4Sphere.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4PVParameterised.hh"
#include "G4SDManager.hh"

#include "G4UserLimits.hh"

#include "G4VisAttributes.hh"
#include "G4Colour.hh"

#include "G4ios.hh"

#include "G4UnitsTable.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
 
SCDetectorConstruction::SCDetectorConstruction()
:solidWorld(0),  logicWorld(0),  physiWorld(0),
 solidTracker(0),logicTracker(0),physiTracker(0), 
 fpMagField(0), fWorldLength(0.),  fTrackerpDz(0.)
{
  fpMagField = new SCMagneticField();
  detectorMessenger = new SCDetectorMessenger(this);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
 
SCDetectorConstruction::~SCDetectorConstruction()
{
  delete fpMagField;
  delete detectorMessenger;             
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
 
G4VPhysicalVolume* SCDetectorConstruction::Construct()
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

  // Print all the materials defined.
  //
  G4cout << G4endl << "The materials defined are : " << G4endl << G4endl;
  G4cout << *(G4Material::GetMaterialTable()) << G4endl;

//--------- Sizes of the principal geometrical components (solids)  ---------
  
  fWorldLength= 10*m ;
   
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
  // Twisted Tracker
  //------------------------------
  
  G4ThreeVector positionTracker = G4ThreeVector(0,0,0);
  G4RotationMatrix* rotD3 = new G4RotationMatrix();

  /*
  fTwistAngle = 60*deg ;    // twisted angle
  fTrackerpDx1 = 15*cm ;
  fTrackerpDx2 = 25*cm ;
  fTrackerpDy1 = 30*cm ;
  fTrackerpDx3 = 15*cm ;
  fTrackerpDx4 = 25*cm ;
  fTrackerpDy2 = 30*cm ;
  fAlph = 0*deg ;
  fTheta = 0*deg ;
  fPhi  = 0*deg ;
  fTrackerpDz = 80*cm ;  // half z-length
  */

  /*

  fTwistAngle = 60*deg ;    // twisted angle
  fTrackerpDx1 = 14*cm ;
  fTrackerpDx2 = 24*cm ;
  fTrackerpDy1 = 16*cm ;
  fTrackerpDx3 = 11*cm ;
  fTrackerpDx4 = 16*cm ;
  fTrackerpDy2 = 8*cm ;
  fAlph = 40*deg ;
  fTheta = 10*deg ;
  fPhi  = 40*deg ;
  fTrackerpDz = 80*cm ;  // half z-length

  */

  fTwistAngle = 20*deg ;
  fTrackerpDz = 80*cm ;
  fTheta = 10*deg ;
  fPhi  =  40*deg ;
  fTrackerpDy1 = 8*cm ;
  fTrackerpDx1 = 11*cm ;
  fTrackerpDx2 = 16*cm ;
  fTrackerpDy2 = 8*cm ;
  fTrackerpDx3 = 11*cm ;
  fTrackerpDx4 = 16*cm ;
  fAlph = -50*deg    ;

  /*
  solidTracker = new G4TwistedTrap("tracker",
				      fTwistAngle,    // twist angle
				      fTrackerpDz,         // half z length
				      fTheta,      // direction between end planes
				      fPhi,        // defined by polar and azimutal angles.
				      fTrackerpDy1,        // half y length at -pDz
				      fTrackerpDx1,        // half x length at -pDz,-pDy
				      fTrackerpDx2,        // half x length at -pDz,+pDy
				      fTrackerpDy2,        // half y length at +pDz
				      fTrackerpDx3,         // half x length at +pDz,-pDy
				      fTrackerpDx4,        // half x length at +pDz,+pDy
				      -fAlph        // tilt angle . attention: inverse definition.
				      ) ;


  */


  solidTracker = new G4Torus("tracker", 0, fTrackerpDy2 ,fTrackerpDx2, 0*deg, 360*deg) ;


  logicTracker = new G4LogicalVolume(solidTracker , Air, "Tracker",0,0,0);  
  physiTracker = new G4PVPlacement(rotD3,         // no rotation
				  positionTracker, // at (x,y,z)
				  logicTracker,    // its logical volume		  
				  "Tracker",       // its name
				  logicWorld,      // its mother  volume
				  false,           // no boolean operations
				  0);              // no particular field 

  G4cout << "Mass of twisted Trapezoid = " << G4BestUnit(solidTracker->GetCubicVolume(),"Volume") << G4endl ;

//--------- Visualization attributes -------------------------------


// the world is transparent
  G4VisAttributes* WorldAtt = new G4VisAttributes(G4Colour(1.,1.,1.,0.));
  WorldAtt->SetVisibility(true);
  logicWorld->SetVisAttributes(WorldAtt);  

  G4VisAttributes* BoxVisAtt= new G4VisAttributes(G4Colour(0.0,.0,1.0,0.6));
  BoxVisAtt->SetVisibility(true);
  logicTracker->SetVisAttributes(BoxVisAtt);
  
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
 
 
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
 
void SCDetectorConstruction::SetMagField(G4double fieldValue)
{
  fpMagField->SetFieldValue(fieldValue);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
