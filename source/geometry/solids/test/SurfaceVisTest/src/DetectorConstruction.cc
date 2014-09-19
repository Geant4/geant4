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
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo...... 
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
 
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include "DetectorConstruction.hh"
#include "DetectorMessenger.hh"

#include "G4Material.hh"

#include "G4Box.hh"
#include "G4Orb.hh"
#include "G4Tubs.hh"
#include "G4CutTubs.hh"
#include "G4Sphere.hh"
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

#include "G4BooleanSolid.hh"
#include "G4DisplacedSolid.hh"
#include "G4UnionSolid.hh"
#include "G4IntersectionSolid.hh"
#include "G4SubtractionSolid.hh"
#include "G4ReflectedSolid.hh"

#include "G4TwoVector.hh"
#include "G4TriangularFacet.hh"
#include "G4QuadrangularFacet.hh"

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
#include "G4GeometryManager.hh"
#include "G4PhysicalVolumeStore.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4SolidStore.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
 
DetectorConstruction::DetectorConstruction()
:solidWorld(0),  logicWorld(0),  physiWorld(0),
 physiTracker(0), fWorldLength(0.)
{
    detectorMessenger = new DetectorMessenger(this);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
 
DetectorConstruction::~DetectorConstruction()
{
    delete detectorMessenger;             
}

////////////////////////////////////////////////////////////////

void DetectorConstruction::SwitchDetector()
{
   G4RunManager::GetRunManager()->DefineWorldVolume(physiWorld);
  
}


G4VPhysicalVolume*
DetectorConstruction::SelectDetector( const G4String& val )
{

   G4GeometryManager::GetInstance()->OpenGeometry();
   G4PhysicalVolumeStore::GetInstance()->Clean();
   G4LogicalVolumeStore::GetInstance()->Clean();
   G4SolidStore::GetInstance()->Clean();
   ///////
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
  G4double fTrackerR1,fTrackerR2,fPhi,fPhiSegment,fTheta,fThetaSegment ;  
  G4double fSemiAxisX,fSemiAxisY,fSemiAxisZ  ;
  G4double fTrackerR,fTrackerpDx1,fTrackerpDy1,fTrackerpDz;
  if (val == "Sphere")
  {
    // only solid sphere is supportet for now
      fTrackerR1 = 0*cm ;
      fTrackerR2 = 12*cm ;
      fPhi = 0*deg ;
      fPhiSegment = 360*deg ;
      fTheta = 0*deg ;
      fThetaSegment = 180*deg ;
      aVolume = new G4Sphere ("aSphere",  fTrackerR1, fTrackerR2 ,fPhi, fPhiSegment, fTheta, fThetaSegment);
  }

  else if (val == "HalfSphere")
  {
    // only solid sphere is supportet for now
      fTrackerR1 = 0*cm ;
      fTrackerR2 = 12*cm ;
      fPhi = 0*deg ;
      fPhiSegment = 180*deg ;
      fTheta = 0*deg ;
      fThetaSegment = 180*deg ;
      aVolume = new G4Sphere ("aHalfSphere",  fTrackerR1, fTrackerR2 ,fPhi, fPhiSegment, fTheta, fThetaSegment);

  }
  else if (val == "HollowSphere")
  {
    // only solid sphere is supportet for now
      fTrackerR1 = 12*cm ;
      fTrackerR2 = 14*cm ;
      fPhi = 0*deg ;
      fPhiSegment = 360*deg ;
      fTheta = 0*deg ;
      fThetaSegment = 180*deg ;
      aVolume = new G4Sphere ("aHollowSphere",  fTrackerR1, fTrackerR2 ,fPhi, fPhiSegment, fTheta, fThetaSegment);

  }
  else if (val == "HalfHollowSphere")
  {
    // only solid sphere is supportet for now
      fTrackerR1 = 12*cm ;
      fTrackerR2 = 14*cm ;
      fPhi = 0*deg ;
      fPhiSegment = 180*deg ;
      fTheta = 0*deg ;
      fThetaSegment = 180*deg ;
      aVolume = new G4Sphere ("aHalfHollowSphere",  fTrackerR1, fTrackerR2 ,fPhi, fPhiSegment, fTheta, fThetaSegment);

  }
  else if (val == "Q1Shell")
  {
    // only solid sphere is supportet for now
      fTrackerR1 = 10*cm ;
      fTrackerR2 = 12*cm ;
      fPhi = 0*deg ;
      fPhiSegment = 90*deg ;
      fTheta = 0*deg ;
      fThetaSegment = 90*deg ;
      aVolume = new G4Sphere ("aQ1Shell",  fTrackerR1, fTrackerR2 ,fPhi, fPhiSegment, fTheta, fThetaSegment);

  }
  else if (val == "Q2Shell")
  {
    // only solid sphere is supportet for now
      fTrackerR1 = 10*cm ;
      fTrackerR2 = 12*cm ;
      fPhi = 90*deg ;
      fPhiSegment = 90*deg ;
      fTheta = 0*deg ;
      fThetaSegment = 90*deg ;
      aVolume = new G4Sphere ("aQ2Shell",  fTrackerR1, fTrackerR2 ,fPhi, fPhiSegment, fTheta, fThetaSegment);

  }  
  else if (val == "Q3Shell")
  {
    // only solid sphere is supportet for now
      fTrackerR1 = 10*cm ;
      fTrackerR2 = 12*cm ;
      fPhi = 180*deg ;
      fPhiSegment = 90*deg ;
      fTheta = 0*deg ;
      fThetaSegment = 90*deg ;
      aVolume = new G4Sphere ("aQ3Shell",  fTrackerR1, fTrackerR2 ,fPhi, fPhiSegment, fTheta, fThetaSegment);

  }
  else if (val == "Q4Shell")
  {
    // only solid sphere is supportet for now
      fTrackerR1 = 10*cm ;
      fTrackerR2 = 12*cm ;
      fPhi = 270*deg ;
      fPhiSegment = 90*deg ;
      fTheta = 0*deg ;
      fThetaSegment = 90*deg ;
      aVolume = new G4Sphere ("aQ4Shell",  fTrackerR1, fTrackerR2 ,fPhi, fPhiSegment, fTheta, fThetaSegment);

  }
  else if (val == "Q5Shell")
  {
    // only solid sphere is supportet for now
      fTrackerR1 = 10*cm ;
      fTrackerR2 = 12*cm ;
      fPhi = 0*deg ;
      fPhiSegment = 90*deg ;
      fTheta = 90*deg ;
      fThetaSegment = 90*deg ;
      aVolume = new G4Sphere ("aQ5Shell",  fTrackerR1, fTrackerR2 ,fPhi, fPhiSegment, fTheta, fThetaSegment);

  }
else if (val == "Q6Shell")
  {
    // only solid sphere is supportet for now
      fTrackerR1 = 10*cm ;
      fTrackerR2 = 12*cm ;
      fPhi = 90*deg ;
      fPhiSegment = 90*deg ;
      fTheta = 90*deg ;
      fThetaSegment = 90*deg ;
      aVolume = new G4Sphere ("aQ5Shell",  fTrackerR1, fTrackerR2 ,fPhi, fPhiSegment, fTheta, fThetaSegment);

  }  
  else if (val == "Q7Shell")
  {
    // only solid sphere is supportet for now
      fTrackerR1 = 10*cm ;
      fTrackerR2 = 12*cm ;
      fPhi = 180*deg ;
      fPhiSegment = 90*deg ;
      fTheta = 90*deg ;
      fThetaSegment = 90*deg ;
      aVolume = new G4Sphere ("aQ7Shell",  fTrackerR1, fTrackerR2 ,fPhi, fPhiSegment, fTheta, fThetaSegment);

  }
  else if (val == "Q8Shell")
  {
    // only solid sphere is supportet for now
      fTrackerR1 = 10*cm ;
      fTrackerR2 = 12*cm ;
      fPhi = 270*deg ;
      fPhiSegment = 90*deg ;
      fTheta = 90*deg ;
      fThetaSegment = 90*deg ;
      aVolume = new G4Sphere ("aQ8Shell",  fTrackerR1, fTrackerR2 ,fPhi, fPhiSegment, fTheta, fThetaSegment);

  }


  else if (val == "Shell")
  {
    // only solid sphere is supportet for now
      fTrackerR1 = 10*cm ;
      fTrackerR2 = 12*cm ;
      fPhi = 30*deg ;
      fPhiSegment = 120*deg ;
      fTheta = 10*deg ;
      fThetaSegment = 100*deg ;
      aVolume = new G4Sphere ("aShell",  fTrackerR1, fTrackerR2 ,fPhi, fPhiSegment, fTheta, fThetaSegment);

  }
  else if (val == "Ring")
  {
    // only solid sphere is supportet for now
      fTrackerR1 = 10*cm ;
      fTrackerR2 = 12*cm ;
      fPhi = 30*deg ;
      fPhiSegment = 120*deg ;
      fTheta = 10*deg ;
      fThetaSegment = 40*deg ;
      aVolume = new G4Sphere ("aRing",  fTrackerR1, fTrackerR2 ,fPhi, fPhiSegment, fTheta, fThetaSegment);

  }
  else if ( val == "Ellipsoid" ) {

    fSemiAxisX = 1*cm ;
    fSemiAxisY = 10*cm ;
    fSemiAxisZ = 100*cm ;

    aVolume = new  G4Ellipsoid("aEllipsoid",
			       fSemiAxisX,
			       fSemiAxisY,
			       fSemiAxisZ,
			       -fSemiAxisZ,fSemiAxisZ);


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
  else if (val == "Cons")
  {        

    fTrackerpDz = 80*cm ;
    fTrackerR1 = 11*cm ;
    fTrackerR2 = 16*cm ;

    aVolume = new G4Cons ( "aCons", 0.8*fTrackerR1 , fTrackerR1, 0.8*fTrackerR2, fTrackerR2, fTrackerpDz, 10*deg, 120.*deg  ); 

  }
  else if (val == "manyCons")
  {        
    aVolume = new G4Cons ( "aCone", 2*cm, 6*cm, 8*cm, 14*cm,
                           10*cm, 10*deg, 300*deg ); 
  //  10*cm, 10*deg, 300*deg ); 
			   //  0., pi);

 
  }
  else if (val == "Tubs")
  {

    // only solid Tube is supported.
    fTrackerpDz = 80*cm ;
    fTrackerR = 11*cm ;

    aVolume = new G4Tubs ( "aTube",0.8*fTrackerR,fTrackerR,fTrackerpDz,40.*deg,100*deg) ;

  }
 else if (val == "CutTubs")
  {

    // only solid Tube is supported.
    fTrackerpDz = 80*cm ;
    fTrackerR = 11*cm ;

    aVolume = new G4CutTubs ( "aTube",0.8*fTrackerR,fTrackerR,fTrackerpDz,40.*deg,100*deg,G4ThreeVector(0,-1,-1),G4ThreeVector(0,0.7,0.7)) ;

  }
  else if (val == "Hype")
  {
    aVolume = new G4Hype ("aHype", 7*cm, 10*cm, 40*deg, 40*deg, 40*cm );
  }
  else if (val == "Torus")
  {

    fPhi = 40*deg ;
    fPhiSegment = 100*deg ;

    fTrackerR1 = 5*cm ;
    fTrackerR2 = 6*cm ;
    fTrackerR  = 20*cm ;

    aVolume = new G4Torus("aTorus", fTrackerR1, fTrackerR2 ,fTrackerR, fPhi, fPhiSegment) ;

  }
  else if (val == "Para")
  {
    aVolume = new G4Para ("aPara", 8*cm, 10*cm, 12*cm, 30*deg, 45*deg, 60*deg);
  }
  else if (val == "Paraboloid")
  {
    aVolume = new G4Paraboloid ("aParaboloid", 8*cm, 1*cm, 12*cm);
  }
  else if (val == "PolyconeGen")
  {
    G4double rr[]={100,120,140,120};
    G4double zz[]={-100,-100,100,100};
    aVolume = new G4Polycone ("aPolycone",0*deg,90*deg, 4, rr,zz);
 }
 else if (val == "PolyconeGenComplex")
  {
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
 }
  else if (val == "Polycone")
  {
    G4double rIn[]={100,120,0,0,120,120};
    G4double rOut[]={120,140,140,120,140,160};
    G4double zPlane[]={-200,-100,0,100,200,300};
    aVolume = new G4Polycone ("aPolycone",0*deg,90*deg, 6,zPlane, rIn,rOut);
  }
 else if (val == "PolyhedraGen")
  {
    G4double rr[]={100,120,140,120};
    G4double zz[]={-100,-100,100,100};
    aVolume = new G4Polyhedra ("aPolycone",0*deg,90*deg,8, 4, rr,zz);
  }
 else if (val == "PolyhedraGenComplex")
  {
    G4double rr[]={ 100, 120,   0,  0,140,120,150,150,150,140, 130};
    G4double zz[]={-150,-100,-100,-50,  0,200,200,100,0,  -50,-100};
    aVolume = new G4Polyhedra ("aPolycone",0*deg,360*deg,8, 11, rr,zz);
  }
 else if (val == "Polyhedra")
  {
    G4double rOut[5]={2100,600,1500,1500,3800.};
    G4double rIn[5]={0,500,0,700,100};
    G4double zPlane[5]={-100,1000,1500,2500,3000.};
    aVolume = new G4Polyhedra ("aPolycone",0*deg,60*deg,6, 5,zPlane, rIn,rOut);
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

  else if ( val == "Tet" ) 
  {
    //   G4ThreeVector anchor = G4ThreeVector(0, 0, 0);
    //G4ThreeVector     p2 = G4ThreeVector(10*cm, 5*cm , 0);
    //G4ThreeVector     p3 = G4ThreeVector(5*cm,10*cm,0);
    //G4ThreeVector     p4 = G4ThreeVector(5*cm,5*cm  ,10*cm);
   
    //aVolume = new G4Tet("aTet",anchor,p2,p3,p4);
      G4ThreeVector pzero(0,0,0);
      G4ThreeVector pnt1(10.*cm,0.*cm,0.*cm),pnt2(5.0*cm,10.*cm,0.*cm), pnt3(5.*cm,5.*cm,10.*cm);
      G4bool  goodTet;
      aVolume= new G4Tet( "aTet", pzero, pnt1, pnt2, pnt3, &goodTet);
  }
  else if ( val == "Trap") 
  {
      fTrackerpDz = 80*cm;
    fTheta = 10*deg ;
    fPhi  =  40*deg ;
    fTrackerpDy1 = 16*cm ;
    fTrackerpDx1 = 24*cm ;
    G4double fTrackerpDx2 = 14*cm ;
    G4double fTrackerpDy2 = 8*cm ;
    G4double fTrackerpDx3 = 16*cm ;
    G4double fTrackerpDx4 = 11*cm ;
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
                        0.5*mm,       // xSemiAxis
                        1*mm,       // ySemiAxis
                        40*mm,      // zheight
                        25*mm) ;    // zTopCut

  }
  else if ( val == "EllipticalTube" ) 
  {
    aVolume = new G4EllipticalTube("aEllipticalTube" ,
				   2*cm,   // xSemiAxis
				   5*cm,   // ySemiAxis
				   35*cm) ;  // zheight

  }
  else if (val == "GenericTrap" ){
   std::vector<G4TwoVector> vertices;
   vertices.push_back(G4TwoVector( -4.5*cm, -4.5*cm));
   vertices.push_back(G4TwoVector( -4.5*cm,  4.5*cm));
   vertices.push_back(G4TwoVector(  4.5*cm,  4.5*cm));
   vertices.push_back(G4TwoVector(  4.5*cm, -4.5*cm));
   vertices.push_back(G4TwoVector( -3.5*cm, -3.5*cm));
   vertices.push_back(G4TwoVector( -3.5*cm,  3.5*cm));
   vertices.push_back(G4TwoVector(  3.5*cm,  3.5*cm));
   vertices.push_back(G4TwoVector(  3.5*cm, -2.5*cm));     
   aVolume = new G4GenericTrap("aGenTrd",14.*cm,vertices);
  }

 else if(val == "TessellatedSolid")
  { 
    G4double targetSize = 10.*cm;
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
   polygon.push_back(G4TwoVector(-3.*cm, -3.0*cm));
   polygon.push_back(G4TwoVector(-3.*cm,  3.0*cm));
   polygon.push_back(G4TwoVector( 3.*cm,  3.0*cm));
   polygon.push_back(G4TwoVector( 3.*cm, -3.0*cm));
   polygon.push_back(G4TwoVector( 1.5*cm, -3.0*cm));
   polygon.push_back(G4TwoVector( 1.5*cm,  1.5*cm));
   polygon.push_back(G4TwoVector(-1.5*cm,  1.5*cm));
   polygon.push_back(G4TwoVector(-1.5*cm, -3.0*cm));
  
   std::vector<G4ExtrudedSolid::ZSection> zsections;
   zsections.push_back(G4ExtrudedSolid::ZSection(-4.0*cm, G4TwoVector(-2.0*cm, 1.0*cm), 1.5));
   zsections.push_back(G4ExtrudedSolid::ZSection( 1.0*cm, G4TwoVector(  0*cm,  0*cm), 0.5));
   zsections.push_back(G4ExtrudedSolid::ZSection( 1.5*cm, G4TwoVector(  0*cm,  0*cm), 0.7));
   zsections.push_back(G4ExtrudedSolid::ZSection( 4.0*cm, G4TwoVector( 2.0*cm, 2.0*cm), 0.9));

   aVolume = new G4ExtrudedSolid("aExtrudedSolid", polygon, zsections);
  }
   else if (val == "TwistedBox")
  {
    aVolume = new G4TwistedBox("aTwistedBox",40*deg,5*cm,10*cm,15*cm);
  }
  else if (val == "TwistedTrd")
  { 
    aVolume = new G4TwistedTrd("aTwistedTrd",5*cm,10*cm,8*cm,15*cm,18*cm,20*deg);
  }
  else if (val == "TwistedTrap")
  {
    aVolume = new G4TwistedTrap("aTwistedTrap",40*deg,5*cm,10*cm,8*cm,15*cm);
  }
  else if ( val == "TwistedTrap2") 
  {
    aVolume = new G4TwistedTrap("aTwistedTrap2",
				   20*deg,    // twist angle
				   80*cm,         // half z length
				   10*deg,      // direction between end planes
				   40*deg,        // defined by polar and azimutal angles.
				   8*cm,        // half y length at -pDz
				   11*cm,        // half x length at -pDz,-pDy
				   16*cm,        // half x length at -pDz,+pDy
				   8*cm,        // half y length at +pDz
				   11*cm,         // half x length at +pDz,-pDy
				   16*cm,        // half x length at +pDz,+pDy
				   -50*deg        // tilt angle at +pDz
				   ) ;
  }
  else if ( val == "TwistedTubs")
  {
    aVolume = new G4TwistedTubs("aTwistedTubs",10.*deg,1*cm,2*cm,4*cm,171.*deg);

  }
  else
    { G4ExceptionDescription desc;
      desc << "DetectorConstruction tried to select " << val << G4endl;
      G4Exception("DetectorConstruction::SelectDetector()", "SVT001",
		  FatalException, desc);
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

  // G4VPhysicalVolume * aVolume_phys1 =
      new G4PVPlacement(0,
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

  G4VisAttributes* BoxVisAtt= new G4VisAttributes(G4Colour(0.0,1.0,0.0));
  BoxVisAtt->SetVisibility(true);
  BoxVisAtt->SetForceSolid(true);
  aVolume_log->SetVisAttributes(BoxVisAtt);
  
//--------- example of User Limits -------------------------------

  
  return physiWorld;

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
 
G4VPhysicalVolume* DetectorConstruction::Construct()
{


  //-------------------Hall ----------------------------------
  
  return SelectDetector ("Polyhedra");  // default is Polyhedra

}

