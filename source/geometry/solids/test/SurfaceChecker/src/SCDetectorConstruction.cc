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
// $Id: SCDetectorConstruction.cc,v 1.3 2005-07-01 12:13:50 link Exp $
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
#include "G4Orb.hh"
#include "G4Tubs.hh"
#include "G4Sphere.hh"
#include "G4Cons.hh"
#include "G4Hype.hh"
#include "G4Para.hh"
#include "G4Torus.hh"
#include "G4Trd.hh"

#include "G4Polycone.hh"

#include "G4TwistedTubs.hh"
#include "G4TwistedBox.hh"
#include "G4TwistedTrd.hh"
#include "G4TwistedTrap.hh"

#include "G4BooleanSolid.hh"
#include "G4DisplacedSolid.hh"
#include "G4UnionSolid.hh"
#include "G4IntersectionSolid.hh"
#include "G4SubtractionSolid.hh"
#include "G4ReflectedSolid.hh"

#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4PVParameterised.hh"
#include "G4SDManager.hh"

#include "G4UserLimits.hh"

#include "G4VisAttributes.hh"
#include "G4Colour.hh"

#include "G4ios.hh"

#include "G4UnitsTable.hh"

#include "G4RunManager.hh"


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
 
SCDetectorConstruction::SCDetectorConstruction()
:solidWorld(0),  logicWorld(0),  physiWorld(0),
 logicTracker(0),physiTracker(0), 
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

////////////////////////////////////////////////////////////////

void SCDetectorConstruction::SwitchDetector()
{
  G4RunManager::GetRunManager()->DefineWorldVolume(physiWorld);
}


G4VPhysicalVolume*
SCDetectorConstruction::SelectDetector( const G4String& val )
{


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
  //  G4cout << G4endl << "The materials defined are : " << G4endl << G4endl;
  //G4cout << *(G4Material::GetMaterialTable()) << G4endl;

  G4Box* b1 = new G4Box ( "b1", 100*cm, 50*cm, 50*cm );
  G4Box* b2 = new G4Box ( "b2", 50*cm, 100*cm, 50*cm );

  fval = val ;

  if (val == "Sphere")
  {
    // only solid sphere is supportet for now
      fTrackerR1 = 10*cm ;
      fTrackerR2 = 12*cm ;
      fPhi = 30*deg ;
      fPhiSegment = 120*deg ;
      fTheta = 10*deg ;
      fThetaSegment = 100*deg ;
      aVolume = new G4Sphere ("aSphere",  fTrackerR1, fTrackerR2 ,fPhi, fPhiSegment, fTheta, fThetaSegment);

  }
  else if (val == "Orb")
  {

    fTrackerR = 11*cm ;
    aVolume = new G4Orb ( "aOrb", fTrackerR );

  }
  else if (val == "Box") 
  {         

    fTrackerpDx1 = 10*cm ;
    fTrackerpDy1 = 10*cm ;
    fTrackerpDz  = 10*cm ;

    aVolume = new G4Box ( "aBox", fTrackerpDx1, fTrackerpDy1, fTrackerpDz );
  }
  else if (val == "Cone")
  {        

    fTrackerpDz = 80*cm ;
    fTrackerR1 = 11*cm ;
    fTrackerR2 = 16*cm ;

    aVolume = new G4Cons ( "aCone", 0 , fTrackerR1, 0, fTrackerR2, fTrackerpDz, 0, 360.*deg  ); 

  }
  else if (val == "manyCons")
  {        
    aVolume = new G4Cons ( "aCone", 2*cm, 6*cm, 8*cm, 14*cm,
                           10*cm, 10*deg, 300*deg ); 
  //  10*cm, 10*deg, 300*deg ); 
			   //  0., pi);

 
  }
  else if (val == "Tube")
  {

    // only solid Tube is supported.
    fTrackerpDz = 80*cm ;
    fTrackerR = 11*cm ;

    aVolume = new G4Tubs ( "aTube",0.,fTrackerR,fTrackerpDz,0.,360*deg) ;

  }
  else if (val == "Hype")
  {
    aVolume = new G4Hype ("aHype", 10*cm, 20*cm, 0*deg, 360*deg, 10*cm );
  }
  else if (val == "Torus")
  {

    fTrackerR = 20*cm ;
    fTrackerR1 = 5*cm ;
    fTrackerR2 = 6*cm ;

    aVolume = new G4Torus("aTorus", fTrackerR1, fTrackerR2 ,fTrackerR, 0*deg, 360*deg) ;

  }
  else if (val == "Para")
  {
    aVolume = new G4Para ("aPara", 8*cm, 10*cm, 12*cm, 30*deg, 45*deg, 60*deg);
  }
  else if (val == "Trd")
  {
    aVolume = new G4Trd ("aTrd", 80*cm, 100*cm, 70*cm, 90*cm, 100*cm);
  }
  else if (val == "b1Ub2") 
  {         
    aVolume = new G4UnionSolid("b1Ub2",b1,b2);
    /*
    G4Box * box1 = new G4Box("Box1",1092.500000,240.103374,92.000000);
    G4Box * box2 = new G4Box("Box2",540.103374,792.500000,92.000000);

    G4double L1 = 1104;

    aVolume =
    new G4UnionSolid("ECShapeBoxes",
                     box1,
                     box2,
                     0,
                     G4ThreeVector(-L1/2.,
                                   L1/2.,
                                   0.));
    */
  }
  else if (val == "b1Ib2") 
  {         
    aVolume = new G4IntersectionSolid("b1Ib2",b1,b2);
  }
  else if (val == "b1Sb2") 
  {         
    aVolume = new G4SubtractionSolid("b1Sb2",b1,b2);
  }
  else if (val == "b1Ib1") 
  {         
    aVolume = new G4IntersectionSolid("b1Ib1",b1,b1);
  }
  else if (val == "b1Ub1") 
  {         
    aVolume = new G4UnionSolid("b1Ub1",b1,b1);
  }
  else if (val == "b1Sb1") 
  {         
    aVolume = new G4SubtractionSolid("b1Sb1",b1,b1);
  }
  else if ( val == "TwistedTubs" )
  {

    fTwistAngle = 20*deg ;
    fTrackerpDz = 80*cm ;
    fTrackerR1  = 5*cm ;
    fTrackerR2  = 10*cm ;
    fPhi        = 50*deg ;

    aVolume = new G4TwistedTubs("aTwistedTubs", fTwistAngle, fTrackerR1 , fTrackerR2, fTrackerpDz, fPhi ) ;
  }

  else if (val == "TwistedBox")
  {


    fTwistAngle = 20*deg ;
    fTrackerpDx1 = 11*cm ;
    fTrackerpDy1 = 8*cm ;
    fTrackerpDz  = 80*cm ;

    aVolume = new G4TwistedBox("aTwistedBox",fTwistAngle,fTrackerpDx1,fTrackerpDy1,fTrackerpDz) ;
  }
  else if (val == "TwistedTrd")
  {

    fTrackerpDx1 = 5*cm ;
    fTrackerpDx2 = 10*cm ;
    fTrackerpDy1 = 8*cm ;
    fTrackerpDy2 = 15*cm ;
    fTrackerpDz  = 80*cm ;
    fTwistAngle = 20*deg ;

    aVolume = new G4TwistedTrd("aTwistedTrd",fTrackerpDx1,fTrackerpDx2,fTrackerpDy1,fTrackerpDy2,fTrackerpDz,fTwistAngle);

  }
  else if (val == "TwistedTrap")     // regular twisted trap
  {

    fTwistAngle = 20*deg ;
    fTrackerpDz = 80*cm ;
    fTrackerpDx1 = 5*cm ;
    fTrackerpDx2 = 7*cm ;
    fTrackerpDy1 = 8*cm ;

    aVolume = new G4TwistedTrap("aTwistedTrap",fTwistAngle,fTrackerpDx1,fTrackerpDx2,fTrackerpDy1,fTrackerpDz);
  }

  else if ( val == "TwistedTrap2") 
  {

    // this is a general twisted trap with equal endcaps, but alph,phi,theta not zero
    // Reason: the surface equation is a special case.

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
    fAlph = 50*deg    ;

    aVolume = new G4TwistedTrap("aTwistedTrap2",
				fTwistAngle,         // twist angle
				fTrackerpDz,         // half z length
				fTheta,              // direction between end planes
				fPhi,                // defined by polar and azimutal angles.
				fTrackerpDy1,        // half y length at -pDz
				fTrackerpDx1,        // half x length at -pDz,-pDy
				fTrackerpDx2,        // half x length at -pDz,+pDy
				fTrackerpDy2,        // half y length at +pDz
				fTrackerpDx3,        // half x length at +pDz,-pDy
				fTrackerpDx4,        // half x length at +pDz,+pDy
				fAlph                // tilt angle at +pDz
				) ;
  }
  else if ( val == "TwistedTrap3") 
  {
    fTwistAngle = 60*deg ; 
    fTrackerpDz = 80*cm;
    fTheta = 10*deg ;
    fPhi  =  40*deg ;
    fTrackerpDy1 = 16*cm ;
    fTrackerpDx1 = 24*cm ;
    fTrackerpDx2 = 14*cm ;
    fTrackerpDy2 = 8*cm ;
    fTrackerpDx3 = 16*cm ;
    fTrackerpDx4 = 11*cm ;
    fAlph = 50*deg    ;

    aVolume = new G4TwistedTrap("aTwistedTrap3",
				fTwistAngle,         // twist angle
				fTrackerpDz,         // half z length
				fTheta,              // direction between end planes
				fPhi,                // defined by polar and azimutal angles.
				fTrackerpDy1,        // half y length at -pDz
				fTrackerpDx1,        // half x length at -pDz,-pDy
				fTrackerpDx2,        // half x length at -pDz,+pDy
				fTrackerpDy2,        // half y length at +pDz
				fTrackerpDx3,        // half x length at +pDz,-pDy
				fTrackerpDx4,        // half x length at +pDz,+pDy
				fAlph                // tilt angle at +pDz
				) ;
  }
  else
  {
    G4Exception("Sc01DetectorConstruction::SelectDetector() - Invalid shape!");
  }

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
				 


  G4LogicalVolume* aVolume_log = new G4LogicalVolume(aVolume, Air, "aVolume_L", 0,0,0);

  G4VPhysicalVolume * aVolume_phys1
    = new G4PVPlacement(0,
			G4ThreeVector(0*cm, 0*cm, 0*cm),
                        aVolume_log, 
			val, 
			logicWorld, 
			false,
			0);



//--------- Visualization attributes -------------------------------


// the world is transparent
  G4VisAttributes* WorldAtt = new G4VisAttributes(G4Colour(1.,1.,1.,0.));
  WorldAtt->SetVisibility(true);
  logicWorld->SetVisAttributes(WorldAtt);  

  G4VisAttributes* BoxVisAtt= new G4VisAttributes(G4Colour(0.0,.0,1.0,0.6));
  BoxVisAtt->SetVisibility(true);
  aVolume_log->SetVisAttributes(BoxVisAtt);
  
//--------- example of User Limits -------------------------------

  
  return physiWorld;

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
 
G4VPhysicalVolume* SCDetectorConstruction::Construct()
{


  //-------------------Hall ----------------------------------
  
  return SelectDetector ("Sphere");  // default is Sphere

}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
 
void SCDetectorConstruction::SetMagField(G4double fieldValue)
{
  fpMagField->SetFieldValue(fieldValue);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
