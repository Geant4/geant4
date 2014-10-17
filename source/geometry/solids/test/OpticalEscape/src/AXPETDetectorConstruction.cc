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
// $Id$
// ------------------------------------------------------------
// Geant4 class implementation file
//
// 03/09/2008, by T.Nikitina
// ------------------------------------------------------------

#include "AXPETDetectorConstruction.hh"
#include "AXPETDetectorMessenger.hh"

#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"

//Solids
#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4CutTubs.hh"
#include "G4Sphere.hh"
#include "G4Orb.hh"
#include "G4Cons.hh"
#include "G4Hype.hh"
#include "G4Para.hh"
#include "G4Paraboloid.hh"
#include "G4Torus.hh"
#include "G4Trd.hh"
#include "G4Ellipsoid.hh"
#include "G4EllipticalTube.hh"
#include "G4EllipticalCone.hh"

#include "G4Trap.hh"
#include "G4Trd.hh"
#include "G4Tet.hh"
#include "G4GenericTrap.hh"
#include "G4TessellatedSolid.hh"
#include "G4ExtrudedSolid.hh"

#include "G4Polycone.hh"
#include "G4Polyhedra.hh"

#include "G4TwistedTubs.hh"
#include "G4TwistedBox.hh"
#include "G4TwistedTrd.hh"
#include "G4TwistedTrap.hh"

#include "G4TwoVector.hh"
#include "G4TriangularFacet.hh"
#include "G4QuadrangularFacet.hh"

#include "G4BooleanSolid.hh"
#include "G4DisplacedSolid.hh"
#include "G4UnionSolid.hh"
#include "G4IntersectionSolid.hh"
#include "G4SubtractionSolid.hh"
#include "G4ReflectedSolid.hh"
#include "G4RotationMatrix.hh"
//
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4PVReplica.hh"
#include "G4ThreeVector.hh"
#include "G4RotationMatrix.hh"

#include "G4RunManager.hh"
#include "G4GeometryManager.hh"
#include "G4PhysicalVolumeStore.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4SolidStore.hh"

#include "G4MaterialPropertiesTable.hh"
#include "G4OpticalSurface.hh"
#include "G4LogicalSkinSurface.hh"
#include "G4LogicalBorderSurface.hh"
#include "AXPETMaterial.hh"

#include "G4VisAttributes.hh"
#include "G4Colour.hh"

#include "globals.hh"

AXPETDetectorConstruction::AXPETDetectorConstruction():
WorldVolume(0),LogWorldVolume(0),PhysWorldVolume(0), 
						      LYSOVolume(0),LogLYSOVolume(0), PhysLYSOVolume(0)//, 
  //InterVolume(0),LogInterVolume(0), PhysInterVolume(0), 
  //DetecVolume(0),LogDetecVolume(0), PhysDetecVolume(0) 
{
 pttoMaterial = new AXPETMaterial();

 worldXsize  = 20. *cm;
 worldYsize  = 20. *cm;
 worldZsize  = 20. *cm;

 lysoIz   =  2.0 *mm;
 lysoDz   =  0.1 *mm;
 lysoR    =  6.0 *mm;

 detectorMessenger = new AXPETDetectorMessenger(this);
 fval="Tubs";

 xRot=0.0;
 yRot=0.0;
 zRot=0.0;
 fAbort=0;
}

AXPETDetectorConstruction::~AXPETDetectorConstruction()
{    delete detectorMessenger;     }



void AXPETDetectorConstruction::SwitchDetector()
{
   
   G4RunManager::GetRunManager()->DefineWorldVolume(PhysWorldVolume);
  
}
G4VSolid*
AXPETDetectorConstruction::SelectDetector( const G4String& val )
{

  G4Box* b1 = new G4Box ( "b1", 1*mm, 0.5*mm, 0.50*mm );
  G4Box* b2 = new G4Box ( "b2", 0.5*mm, 1*mm, 0.5*mm );

  fval = val ;
  G4double fTrackerR1,fTrackerR2,fPhi,fPhiSegment,fTheta,fThetaSegment ;  
  G4double fSemiAxisX,fSemiAxisY,fSemiAxisZ  ;
  G4double fTrackerR,fTrackerpDx1,fTrackerpDy1,fTrackerpDz;
  if (val == "Sphere")
  {
    // only solid sphere is supportet for now
      fTrackerR1 = 0*mm ;
      fTrackerR2 = (lysoR-0.01*mm)/2. ;
      fPhi = 0*deg ;
      fPhiSegment = 360*deg ;
      fTheta = 0*deg ;
      fThetaSegment = 180*deg ;
      aVolume = new G4Sphere ("aSphere",  fTrackerR1, fTrackerR2 ,fPhi, fPhiSegment, fTheta, fThetaSegment);
  }

  else if (val == "HalfSphere")
  {
    // only solid sphere is supportet for now
      fTrackerR1 = 0. ;
      fTrackerR2 = (lysoR-0.01*mm)/2. ;
      fPhi = 0*deg ;
      fPhiSegment = 180*deg ;
      fTheta = 0*deg ;
      fThetaSegment = 180*deg ;
      aVolume = new G4Sphere ("aHalfSphere",  fTrackerR1, fTrackerR2 ,fPhi, fPhiSegment, fTheta, fThetaSegment);

  }
  else if (val == "HollowSphere")
  {
    // only solid sphere is supportet for now
      fTrackerR1 = (lysoR-0.01*mm)/4.;
      fTrackerR2 = (lysoR-0.01*mm)/2. ;
      fPhi = 0*deg ;
      fPhiSegment = 360*deg ;
      fTheta = 0*deg ;
      fThetaSegment = 180*deg ;
      aVolume = new G4Sphere ("aHollowSphere",  fTrackerR1, fTrackerR2 ,fPhi, fPhiSegment, fTheta, fThetaSegment);

  }
  else if (val == "HalfHollowSphere")
  {
    // only solid sphere is supportet for now
      fTrackerR1 =(lysoR-0.01*mm)/4.  ;
      fTrackerR2 =(lysoR-0.01*mm)/2. ;
      fPhi = 0*deg ;
      fPhiSegment = 180*deg ;
      fTheta = 0*deg ;
      fThetaSegment = 180*deg ;
      aVolume = new G4Sphere ("aHalfHollowSphere",  fTrackerR1, fTrackerR2 ,fPhi, fPhiSegment, fTheta, fThetaSegment);

  }
  else if (val == "Q1Shell")
  {
    // only solid sphere is supportet for now
      fTrackerR1 =(lysoR-0.01*mm)/4.  ;
      fTrackerR2 =(lysoR-0.01*mm)/2. ;
      fPhi = 0*deg ;
      fPhiSegment = 90*deg ;
      fTheta = 0*deg ;
      fThetaSegment = 90*deg ;
      aVolume = new G4Sphere ("aQ1Shell",  fTrackerR1, fTrackerR2 ,fPhi, fPhiSegment, fTheta, fThetaSegment);

  }
  else if (val == "Q2Shell")
  {
    // only solid sphere is supportet for now
      fTrackerR1 =(lysoR-0.01*mm)/4.  ;
      fTrackerR2 =(lysoR-0.01*mm)/2. ;
      fPhi = 90*deg ;
      fPhiSegment = 90*deg ;
      fTheta = 0*deg ;
      fThetaSegment = 90*deg ;
      aVolume = new G4Sphere ("aQ2Shell",  fTrackerR1, fTrackerR2 ,fPhi, fPhiSegment, fTheta, fThetaSegment);

  }  
  else if (val == "Q3Shell")
  {
    // only solid sphere is supportet for now
      fTrackerR1 =(lysoR-0.01*mm)/4.  ;
      fTrackerR2 =(lysoR-0.01*mm)/2. ;
      fPhi = 180*deg ;
      fPhiSegment = 90*deg ;
      fTheta = 0*deg ;
      fThetaSegment = 90*deg ;
      aVolume = new G4Sphere ("aQ3Shell",  fTrackerR1, fTrackerR2 ,fPhi, fPhiSegment, fTheta, fThetaSegment);

  }
  else if (val == "Q4Shell")
  {
    // only solid sphere is supportet for now
      fTrackerR1 =(lysoR-0.01*mm)/4.  ;
      fTrackerR2 =(lysoR-0.01*mm)/2. ;
      fPhi = 270*deg ;
      fPhiSegment = 90*deg ;
      fTheta = 0*deg ;
      fThetaSegment = 90*deg ;
      aVolume = new G4Sphere ("aQ4Shell",  fTrackerR1, fTrackerR2 ,fPhi, fPhiSegment, fTheta, fThetaSegment);

  }
  else if (val == "Q5Shell")
  {
    // only solid sphere is supportet for now
       fTrackerR1 =(lysoR-0.01*mm)/4.  ;
      fTrackerR2 =(lysoR-0.01*mm)/2. ;
      fPhi = 0*deg ;
      fPhiSegment = 90*deg ;
      fTheta = 90*deg ;
      fThetaSegment = 90*deg ;
      aVolume = new G4Sphere ("aQ5Shell",  fTrackerR1, fTrackerR2 ,fPhi, fPhiSegment, fTheta, fThetaSegment);

  }
else if (val == "Q6Shell")
  {
    // only solid sphere is supportet for now
      fTrackerR1 =(lysoR-0.01*mm)/4.  ;
      fTrackerR2 =(lysoR-0.01*mm)/2. ;
      fPhi = 90*deg ;
      fPhiSegment = 90*deg ;
      fTheta = 90*deg ;
      fThetaSegment = 90*deg ;
      aVolume = new G4Sphere ("aQ5Shell",  fTrackerR1, fTrackerR2 ,fPhi, fPhiSegment, fTheta, fThetaSegment);

  }  
  else if (val == "Q7Shell")
  {
    // only solid sphere is supportet for now
       fTrackerR1 =(lysoR-0.01*mm)/4.  ;
      fTrackerR2 =(lysoR-0.01*mm)/2. ;
      fPhi = 180*deg ;
      fPhiSegment = 90*deg ;
      fTheta = 90*deg ;
      fThetaSegment = 90*deg ;
      aVolume = new G4Sphere ("aQ7Shell",  fTrackerR1, fTrackerR2 ,fPhi, fPhiSegment, fTheta, fThetaSegment);

  }
  else if (val == "Q8Shell")
  {
    // only solid sphere is supportet for now
      fTrackerR1 =(lysoR-0.01*mm)/4.  ;
      fTrackerR2 =(lysoR-0.01*mm)/2. ;
      fPhi = 270*deg ;
      fPhiSegment = 90*deg ;
      fTheta = 90*deg ;
      fThetaSegment = 90*deg ;
      aVolume = new G4Sphere ("aQ8Shell",  fTrackerR1, fTrackerR2 ,fPhi, fPhiSegment, fTheta, fThetaSegment);

  }


  else if (val == "Shell")
  {
   
 // Test for Problem Report with Optical Photons and Sphere or Ellipsoid
     G4Sphere *newSOut = new G4Sphere("Crystal",0.,lysoR-0.01/mm,0,twopi,0,pi/2);
     G4Sphere *newSIn = new G4Sphere("Crystal",0.,lysoR-4.01/mm,0,twopi,0,twopi);

     // Box was added in user application in order to avoid problem with Sphere with theta=90 deg
     // Not need anymore, this bug in Sphere was repared in 4.9.2.  
     // G4Box* farMirrorBox2 = new G4Box("nearMirrorBox", lysoR-0.01/mm,
     // lysoR-0.01/mm,  lysoR-0.01/mm);
     // aVolume = new G4IntersectionSolid("GlassMirror", farMirrorBox2,shell,0,G4ThreeVector(0,0,-lysoR));
   aVolume = new G4SubtractionSolid("OutsideminInside", newSOut, newSIn);     

  }
  else if (val == "Ring")
  {
    // only solid sphere is supportet for now
      fTrackerR1 =(lysoR-0.01*mm)/4.  ;
      fTrackerR2 =(lysoR-0.01*mm)/2. ;
      fPhi = 30*deg ;
      fPhiSegment = 120*deg ;
      fTheta = 10*deg ;
      fThetaSegment = 40*deg ;
      aVolume = new G4Sphere ("aRing",  fTrackerR1, fTrackerR2 ,fPhi, fPhiSegment, fTheta, fThetaSegment);

  }
  else if ( val == "Ellipsoid" ) {

    fSemiAxisX = (lysoR-0.01*mm)/4. ;
    fSemiAxisY = (lysoR-0.01*mm)/2. ;
    fSemiAxisZ = (lysoR-0.01*mm)/2. ;

    aVolume = new  G4Ellipsoid("aEllipsoid",
			       fSemiAxisX,
			       fSemiAxisY,
			       fSemiAxisZ,
			       -fSemiAxisX,fSemiAxisX);


  }
  else if (val == "Orb")
  {

    fTrackerR = (lysoR-0.01*mm)/2. ;
    aVolume = new G4Orb ( "aOrb", fTrackerR );

  }
  else if (val == "Box") 
  {         

    fTrackerpDx1 = (lysoR-0.01*mm)/2. ;
    fTrackerpDy1 = (lysoR-0.01*mm)/2. ;
    fTrackerpDz  = (lysoR-0.01*mm)/2. ;

    aVolume = new G4Box ( "aBox", fTrackerpDx1, fTrackerpDy1, fTrackerpDz );
  }
  else if (val == "Cons")
  {        

    fTrackerpDz =(lysoR-0.01*mm)/2.  ;
    fTrackerR1 =(lysoR-0.01*mm)/4.  ;
    fTrackerR2 =(lysoR-0.01*mm)/2.  ;

    aVolume = new G4Cons ( "aCons", 0.8*fTrackerR1 , fTrackerR1, 0.8*fTrackerR2, fTrackerR2, fTrackerpDz, 10*deg, 220.*deg  ); 

  }
  else if (val == "manyCons")
  {        
    aVolume = new G4Cons ( "aCone", 2*mm, 3*mm, 2*mm, 4*mm,
                           3.*mm, 10*deg, 300*deg ); 
  //  10*cm, 10*deg, 300*deg ); 
			   //  0., pi);

 
  }
  else if (val == "Tubs")
  {

    // only solid Tube is supported.
  
    aVolume = new G4Tubs ( "aTube",0.1,5.9*mm,3.*mm,0.,180*deg) ;

  }  else if (val == "CutTubs")
  {

     
    aVolume = new G4CutTubs ( "aTube",0.1,5.9*mm,3.*mm,0.,180*deg,G4ThreeVector(0,0,-1),G4ThreeVector(0,0.3,0.95)) ;

  }
  else if (val == "Hype")
  {
    aVolume = new G4Hype ("aHype", 0.7, 1, 40*deg, 40*deg, 4. );
  }
  else if (val == "Torus")
  {

    fPhi = 40*deg ;
    fPhiSegment = 200*deg ;

    fTrackerR1 = 1.*mm ;
    fTrackerR2 = 2.*mm ;
    fTrackerR  = (lysoR-0.01*mm) ;

    aVolume = new G4Torus("aTorus", fTrackerR1, fTrackerR2 ,fTrackerR, fPhi, fPhiSegment) ;

  }
  else if (val == "Para")
  {
    aVolume = new G4Para ("aPara", 0.8, 1.0, 1.2, 30*deg, 45*deg, 60*deg);
  }
  else if (val == "Paraboloid")
  {
    aVolume = new G4Paraboloid ("aParaboloid", 1., 0.5, 1.2);
  }
  else if (val == "PolyconeGen")
  {
    G4double rr[]={1,1.20,1.40,1.20};
    G4double zz[]={-1,-1,1,1};
    aVolume = new G4Polycone ("aPolycone",0*deg,90*deg, 4, rr,zz);
 }
 else if (val == "PolyconeGenComplex")
   { //has to be apdated
    G4double        zPlanes[]       =       {-11640.*mm,    -11640.*mm,     -11570.722*mm,  -9766.005*mm,  
                                                -4487.828*mm,   4487.828*mm,    9766.005*mm,    11570.722*mm,  
                                                11640.*mm,      11640.*mm,      11580.*mm,      9750.*mm,
					    4500.*mm,       -4500.*mm,      -9750.*mm,      -11580.*mm,
					     -11640.*mm};
  
        G4double        radius[]        =       {250.*mm,       305.*mm,        305.*mm,        2770.461*mm,
                                                4932.*mm,       4932.*mm,       2770.461*mm,    305.*mm,
                                                305.*mm,        250.*mm,        250.*mm,        2750.*mm,
						 4900.*mm,       4900.*mm,       2750.*mm,       250.*mm,
						   250.*mm};
	 aVolume = new G4Polycone ("aPolycone",0*deg,360*deg, 16, radius,zPlanes);
    G4VSolid*       SolidPumpPort   =       new     G4Tubs("PumpPort temp tube", 0 * mm, 1750 / 2 * mm, (11161 - 4500) / 2 * mm-1*mm, 0, 360 * deg);
           for (int i=0;i<2;i++)
        {
	       G4RotationMatrix*       Rot     =       new G4RotationMatrix();
	       Rot->rotateZ(45*(-i-1) * deg );
                G4double        X       =       std::cos(45*(i+1)*deg)*-(3050+4750)/2*mm;
                G4double        Y       =       std::sin(45*(i+1)*deg)*-(3050+4750)/2*mm;
                aVolume        =       new G4SubtractionSolid("Main Spectrometer", aVolume, SolidPumpPort,
            Rot, G4ThreeVector(X, Y,(11161+4500)/2*mm - 100 * mm));
	 }



 }
  else if (val == "Polycone")
  {
    G4double rIn[]={1.00,1.20,0,0,1.20,1.20};
    G4double rOut[]={1.20,1.40,1.40,1.20,1.40,1.60};
    G4double zPlane[]={-2.00,-1.00,0,1.00,2.00,3.00};
    aVolume = new G4Polycone ("aPolycone",0*deg,260*deg, 6,zPlane, rIn,rOut);
  }
 else if (val == "PolyhedraGen")
  {
    G4double rr[]={1.00,1.20,1.40,1.20};
    G4double zz[]={-1.00,-1.00,1.00,1.00};
    aVolume = new G4Polyhedra ("aPolycone",0*deg,90*deg,8, 4, rr,zz);
  }
 else if (val == "PolyhedraGenComplex")
  {
    G4double rr[]={ 1.00, 1.20,   0,  0,1.40,1.20,1.50,1.50,1.50,1.40, 1.30};
    G4double zz[]={-1.50,-1.00,-1.00,-0.50,  0,2.00,2.00,1.00,0,  -0.50,-1.00};
    aVolume = new G4Polyhedra ("aPolycone",0*deg,260*deg,8, 11, rr,zz);
  }
 else if (val == "Polyhedra")
  {
    G4double rOut[5]={2.100,0.600,1.500,1.500,3.800};
    G4double rIn[5]={0,0.500,0,0.700,0.100};
    G4double zPlane[5]={-0.100,1.000,1.500,2.500,3.000};
    aVolume = new G4Polyhedra ("aPolycone",0*deg,260*deg,6, 5,zPlane, rIn,rOut);
  }
 
  else if (val == "Trd")
  {
    aVolume = new G4Trd ("aTrd", 0.80, 1.00, 0.70,0.90, 1.00);
  }
  else if (val == "b1Ub2") 
  {         
    aVolume = new G4UnionSolid("b1Ub2",b1,b2);
  
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

  else if ( val == "Tet" ) 
  {
    G4ThreeVector anchor = G4ThreeVector(  0,    0, 0);
    G4ThreeVector     p2 = G4ThreeVector(1.0,  0.5, 0);
    G4ThreeVector     p3 = G4ThreeVector(0.5,  1.0, 0);
    G4ThreeVector     p4 = G4ThreeVector(0.5,  0.5, 1.0);
   
    aVolume = new G4Tet("aTet",anchor,p2,p3,p4);
  }
  else if ( val == "Trap") 
  {
      fTrackerpDz = 0.80;
    fTheta = 10*deg ;
    fPhi  =  40*deg ;
    fTrackerpDy1 = 1.6 ;
    fTrackerpDx1 = 2.4 ;
    G4double fTrackerpDx2 = 1.4 ;
    G4double fTrackerpDy2 = 0.8 ;
    G4double fTrackerpDx3 = 1.6 ;
    G4double fTrackerpDx4 = 1.1 ;
    G4double fAlph = 50*deg    ;

    aVolume = new G4Trap("aTrap",
				fTrackerpDz,         // half z length
				fTheta,              // direction between end planes
				fPhi,                // defined by polar and azimutal angles.
				fTrackerpDy1,        // half y length at -pDz
				fTrackerpDx1,        // half x length at -pDz,-pDy
				fTrackerpDx2,        // half x length at -pDz,+pDy
				fAlph,                // tilt angle at +pDz
				fTrackerpDy2,        // half y length at +pDz
				fTrackerpDx3,        // half x length at +pDz,-pDy
				fTrackerpDx4,        // half x length at +pDz,+pDy
				fAlph                // tilt angle at +pDz
				) ;
    
 
  }
  else if ( val == "EllipticalCone" ) 
  {
    aVolume = new G4EllipticalCone("aEllipticalCone",
                        0.5,       // xSemiAxis
                        1,       // ySemiAxis
                        4.0,      // zheight
                        2.5) ;    // zTopCut

  }
  else if ( val == "EllipticalTube" ) 
  {
    aVolume = new G4EllipticalTube("aEllipticalTube" ,
				   2,   // xSemiAxis
				   4,   // ySemiAxis
				   3.5) ;  // zheight

  }

 else if (val == "TwistedBox")
  {
    aVolume = new G4TwistedBox("aTwistedBox",40*deg,0.5,1.0,1.5);
  }
  else if (val == "TwistedTrd")
  { 
    aVolume = new G4TwistedTrd("aTwistedTrd",0.5,1.0,0.8,1.5,1.8,20*deg);
  }
  else if (val == "TwistedTrap")
  {
    aVolume = new G4TwistedTrap("aTwistedTrap",40*deg,0.5,1.0,0.8,1.5);
  }
  else if ( val == "TwistedTrap2") 
  {
    aVolume = new G4TwistedTrap("aTwistedTrap2",
				   20*deg,    // twist angle
				   0.80,         // half z length
				   10*deg,      // direction between end planes
				   40*deg,        // defined by polar and azimutal angles.
				   0.8,        // half y length at -pDz
				   1.1,        // half x length at -pDz,-pDy
				   1.6,        // half x length at -pDz,+pDy
				   0.8,        // half y length at +pDz
				   1.1,         // half x length at +pDz,-pDy
				   1.6,        // half x length at +pDz,+pDy
				   -50*deg        // tilt angle at +pDz
				   ) ;
  }
  else if ( val == "TwistedTubs")
  {
    aVolume = new G4TwistedTubs("aTwistedTubs",10.*deg,1.0,2.,4.,171.*deg);

  }
 else if (val == "GenericTrap" ){
   std::vector<G4TwoVector> vertices;
   vertices.push_back(G4TwoVector( -4.5, -4.5));
   vertices.push_back(G4TwoVector( -4.5,  4.5));
   vertices.push_back(G4TwoVector(  4.5,  4.5));
   vertices.push_back(G4TwoVector(  4.5, -4.5));
   vertices.push_back(G4TwoVector( -3.5, -3.5));
   vertices.push_back(G4TwoVector( -3.5,  3.5));
   vertices.push_back(G4TwoVector(  3.5,  3.5));
   vertices.push_back(G4TwoVector(  3.5, -2.5));     
   aVolume = new G4GenericTrap("aGenTrd",4.,vertices);
  }
else if(val == "TessellatedSolid")
  { 
    G4double targetSize = 2.;
    G4TessellatedSolid* aVolume1 = new G4TessellatedSolid("aTessellatedSolid");
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

    aVolume1->AddFacet((G4VFacet*) facet1);
    aVolume1->AddFacet((G4VFacet*) facet2);
    aVolume1->AddFacet((G4VFacet*) facet3);
    aVolume1->AddFacet((G4VFacet*) facet4);
    aVolume1->AddFacet((G4VFacet*) facet5);
  
    aVolume1->SetSolidClosed(true);

    aVolume = aVolume1;

  }
  else if (val == "ExtrudedSolid")
  {
   std::vector<G4TwoVector> polygon;
   polygon.push_back(G4TwoVector(-3., -3.0));
   polygon.push_back(G4TwoVector(-3.,  3.0));
   polygon.push_back(G4TwoVector( 3.,  3.0));
   polygon.push_back(G4TwoVector( 3., -3.0));
   polygon.push_back(G4TwoVector( 1.5, -3.0));
   polygon.push_back(G4TwoVector( 1.5,  1.5));
   polygon.push_back(G4TwoVector(-1.5,  1.5));
   polygon.push_back(G4TwoVector(-1.5, -3.0));
  
   std::vector<G4ExtrudedSolid::ZSection> zsections;
   zsections.push_back(G4ExtrudedSolid::ZSection(-4.0, G4TwoVector(-2.0, 1.0), 1.5));
   zsections.push_back(G4ExtrudedSolid::ZSection( 1.0, G4TwoVector(  0,  0), 0.5));
   zsections.push_back(G4ExtrudedSolid::ZSection( 1.5, G4TwoVector(  0,  0), 0.7));
   zsections.push_back(G4ExtrudedSolid::ZSection( 4.0, G4TwoVector( 2.0, 2.0), 0.9));

   aVolume = new G4ExtrudedSolid("aExtrudedSolid", polygon, zsections);
  }


  else
    { G4ExceptionDescription desc;
      desc << "DetectorConstruction tried to select " << val << G4endl;
      G4Exception("DetectorConstruction::SelectDetector()", "AXPET001",
		  FatalException, desc);
    }

  return aVolume;

}
G4VPhysicalVolume* AXPETDetectorConstruction::Construct()
{
 //***********************************************
  // Clear all Stores
  //***********************************************

   G4GeometryManager::GetInstance()->OpenGeometry();
   G4PhysicalVolumeStore::GetInstance()->Clean();
   G4LogicalVolumeStore::GetInstance()->Clean();
   G4SolidStore::GetInstance()->Clean();

 pttoMaterial->DefineMaterials();
 
 air = pttoMaterial->GetMat("AIR");
 lyso= pttoMaterial->GetMat("LYSO");
 vacuum= pttoMaterial->GetMat("Galactic");


  //***********************************************
  // Defining the world volume
  //***********************************************

 WorldVolume = new G4Box("World",worldXsize,worldYsize,worldZsize);

 LogWorldVolume = new G4LogicalVolume(WorldVolume, 
                                           vacuum, 
                                           "logicWorld", 
					    0,0,0);

 PhysWorldVolume = new G4PVPlacement(0, G4ThreeVector(),
					 "physicalWorld", 
					  LogWorldVolume, 
					  0,false,0);


  //***********************************************
  // Defining the lyso-Solid volume
  //***********************************************

 
 LYSOVolume=SelectDetector(fval);
 LogLYSOVolume = new G4LogicalVolume(LYSOVolume, 
                                        lyso, 
                                        "logicLYSO", 
					 0,0,0);

 G4RotationMatrix *mat=new G4RotationMatrix();
 if(xRot!=0.)mat->rotateX(xRot*deg);
 if(yRot!=0.)mat->rotateY(xRot*deg);
 if(zRot!=0.)mat->rotateZ(xRot*deg);
 // mat->rotateY(180.*deg);
 // mat->rotateZ(10.*deg);
 G4ThreeVector position=G4ThreeVector(1.0,1.0,1.0);

 PhysLYSOVolume = new G4PVPlacement(mat,G4ThreeVector(0.0,0.0,0.0),
  					LogLYSOVolume,
					"physicalLYSO", 
				      LogWorldVolume, 
				    	 0,false,0);
					 
  G4cout << "The LYSO crystal is built as "<<fval << G4endl;




  //************************************************************************
  //Al coated surface properties
  //************************************************************************

  const G4int nEntries = 11;
  const G4double refl = 0.85;
  
  G4double PhotonEnergy[nEntries] = {2.478*eV, 2.53*eV, 2.58*eV, 2.636*eV, 2.69*eV, 2.75*eV, 2.816*eV, 2.88*eV, 2.95*eV, 3.022*eV, 3.097*eV};

  G4double Reflectivity[nEntries]   = {refl, refl,refl, refl, refl, refl, refl, refl, refl, refl, refl};
  G4double Efficiency[nEntries]     = {0.0, 0.0,0.0, 0.0,0.0, 0.0,0.0, 0.0,0.0, 0.0,0.0}; 
  G4double specularlobe[nEntries]   = {0.0, 0.0,0.0, 0.0,0.0, 0.0,0.0, 0.0,0.0, 0.0,0.0}; 
  G4double specularspike[nEntries]  = {0.9, 0.9,0.9, 0.9,0.9, 0.9,0.9, 0.9,0.9, 0.9,0.9}; 
  G4double backscatter[nEntries]    = {0.0, 0.0,0.0, 0.0,0.0, 0.0,0.0, 0.0,0.0, 0.0,0.0}; 
  G4double rindex[nEntries]         = {1.39, 1.39,1.39, 1.39,1.39,1.39,1.39, 1.39,1.39, 1.39,1.39}; 


  G4OpticalSurface* alcoated_opsurf=  new G4OpticalSurface("alcoated_opsurf");
   new G4LogicalBorderSurface("alcoated_surf",PhysLYSOVolume,PhysWorldVolume, alcoated_opsurf);
 

  G4MaterialPropertiesTable* alcoated_mt = new G4MaterialPropertiesTable();
  alcoated_mt->AddProperty("RINDEX",  PhotonEnergy,rindex,  nEntries );
  alcoated_mt->AddProperty("EFFICIENCY",  PhotonEnergy,Efficiency,  nEntries );
  alcoated_mt->AddProperty("REFLECTIVITY",PhotonEnergy,Reflectivity,nEntries);
  alcoated_mt->AddProperty("SPECULARLOBECONSTANT",PhotonEnergy,specularlobe,nEntries);
  alcoated_mt->AddProperty("SPECULARSPIKECONSTANT",PhotonEnergy,specularspike,nEntries);
  alcoated_mt->AddProperty("BACKSCATTERCONSTANT",PhotonEnergy,backscatter,nEntries);


  alcoated_opsurf->SetType(dielectric_metal);
  alcoated_opsurf->SetModel(unified);
  alcoated_opsurf->SetFinish(groundbackpainted/*polished*/);
  

  alcoated_opsurf->SetMaterialPropertiesTable(alcoated_mt);


  G4VisAttributes* WhiteVisAtt= new G4VisAttributes(G4Colour(1.0,1.0,1.0));//red
  WhiteVisAtt->SetVisibility(true);
  WhiteVisAtt->SetForceSolid(true);

  G4VisAttributes* LysoVisAtt= new G4VisAttributes(G4Colour(0.0,0.0,1.0));//red
  LysoVisAtt->SetVisibility(true);
  LysoVisAtt->SetForceSolid(true);

  G4VisAttributes* DeteVisAtt= new G4VisAttributes(G4Colour(1.0,1.0,0.0));//yellow
  DeteVisAtt->SetVisibility(true);
  DeteVisAtt->SetForceSolid(true);
  
  LogLYSOVolume->SetVisAttributes(LysoVisAtt);
  LogWorldVolume->SetVisAttributes(G4VisAttributes::GetInvisible());

  return PhysWorldVolume;
 
 }

