//
//
//    ******************************************
//    *                                        *
//    *     ThyroidDetectorConstruction.cc     *
//    *                                        *
//    ******************************************

#include "ThyroidROGeometry.hh"
#include "ThyroidSD.hh"
#include "ThyroidDetectorConstruction.hh"

#include "G4CSGSolid.hh"
#include "G4Sphere.hh"
#include "G4MaterialPropertyVector.hh"
#include "G4SDManager.hh"
#include "G4SubtractionSolid.hh"
#include "G4RunManager.hh"
#include "G4MaterialPropertiesTable.hh"
#include "G4Material.hh"
#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4UnionSolid.hh"
#include "G4LogicalVolume.hh"
#include "G4ThreeVector.hh"
#include "G4PVPlacement.hh"
#include "globals.hh"
#include "G4MaterialTable.hh"
#include "G4Element.hh"
#include "G4ElementTable.hh"
#include "G4PVParameterised.hh"
#include "G4Transform3D.hh"
#include "G4RotationMatrix.hh"
#include "G4FieldManager.hh"
#include "G4TransportationManager.hh"
#include "G4SDManager.hh"
#include "G4VisAttributes.hh"
#include "G4Colour.hh"
#include "G4EllipticalTube.hh"
#include "G4RotationMatrix.hh"

//....

ThyroidDetectorConstruction::ThyroidDetectorConstruction()
{ 

}

//....

ThyroidDetectorConstruction::~ThyroidDetectorConstruction()
{
}

//....

G4VPhysicalVolume* ThyroidDetectorConstruction::Construct()
{

  //definizione colori 
  G4Colour  white   (1.0, 1.0, 1.0) ;
  G4Colour  grey    (0.5, 0.5, 0.5) ;
  G4Colour  lgrey   (.75, .75, .75) ;
  G4Colour  red     (1.0, 0.0, 0.0) ;
  G4Colour  blue    (0.0, 0.0, 1.0) ;
  G4Colour  cyan    (0.0, 1.0, 1.0) ;
  G4Colour  magenta (1.0, 0.0, 1.0) ; 
  G4Colour  yellow  (1.0, 1.0, 0.0) ;
  G4Colour  lblue   (0.0, 0.0, .75);


// Define required materials

 G4double A;    // atomic mass
 G4double Z;    // atomic number
 G4double d;   // density
 G4int iz;     // iz = number of protons in an isotope
 G4int  n;    // n = number of nucleons in an isotope
 G4int ncomponents;
 G4int natoms;
 G4double abundance;
 G4double fractionmass;


 // General elements
 
 A = 1.01*g/mole;
 Z = 1;
 G4Element* elH = new G4Element("Hydrogen","H",Z,A);
  
 A = 14.01*g/mole;
 Z = 7;
 G4Element* elN = new G4Element("Nitrogen","N",Z,A);

 A = 16.00*g/mole;
 Z = 8;
 G4Element* elO = new G4Element("Oxygen","O",Z,A);

 A = 12.01*g/mole;
 Z = 6;
 G4Element* elC = new G4Element("Carbon","C",Z,A);

 A = 23.00*g/mole;
 Z = 11;
 G4Element* elNa = new G4Element("Sodium","Na",Z,A);

 A = 24.30*g/mole;
 Z = 12;
 G4Element* elMg = new G4Element("Magnesium","Mg",Z,A);

 A = 31.00*g/mole;
 Z = 15;
 G4Element* elP = new G4Element("Phosph.","P",Z,A);

 A = 32.00*g/mole;
 Z = 16;
 G4Element* elS = new G4Element("Sulfur","S",Z,A);

 A = 35.45*g/mole;
 Z = 17;
 G4Element* elCl = new G4Element("Chlorine","Cl",Z,A);

 A = 39.10*g/mole;
 Z = 19;
 G4Element* elK = new G4Element("Potassium","K",Z,A);

 A = 40.10*g/mole;
 Z = 20;
 G4Element* elCa = new G4Element("Calcium","Ca",Z,A);

 A = 55.84*g/mole;
 Z = 26;
 G4Element* elFe = new G4Element("Iron","Fe",Z,A);

 A = 65.40*g/mole;
 Z = 30;
 G4Element* elZn = new G4Element("Zinc","Zn",Z,A);
   
 // Air material
 
 d = 1.*mg/cm3;
 G4Material* matAir = new G4Material("Air",d, 2); 
 matAir->AddElement(elN,0.7);
 matAir->AddElement(elO,0.3);

 // SoftTissue

 d = 1.00*g/cm3;
 G4Material* matSoftTissue = new G4Material("SoftTissue",d, ncomponents=13);
 matSoftTissue->AddElement(elH,0.104472);
 matSoftTissue->AddElement(elC,0.232190);
 matSoftTissue->AddElement(elN,0.024880);
 matSoftTissue->AddElement(elO,0.630238);
 matSoftTissue->AddElement(elNa,0.001130);
 matSoftTissue->AddElement(elMg,0.000130);
 matSoftTissue->AddElement(elP,0.001330);
 matSoftTissue->AddElement(elS,0.001990);
 matSoftTissue->AddElement(elCl,0.001340);
 matSoftTissue->AddElement(elK,0.001990);
 matSoftTissue->AddElement(elCa,0.000230);
 matSoftTissue->AddElement(elFe,0.000050);
 matSoftTissue->AddElement(elZn,0.000030);
 
 // Water

 d = 1.000*g/cm3;
 G4Material* matH2O = new G4Material("Water",d,2);
 matH2O->AddElement(elH,2);
 matH2O->AddElement(elO,1);

 // Iodine
 d = 4.93*g/cm3;
 A = 126.90*g/mole;
 Z = 53;
 G4Material* matI131 = new G4Material("I",Z,A,d);

 // Bone

 d=1.85*g/cm3;
 G4Material* bone = new G4Material("bone",d,8);
 bone->AddElement(elH,0.063984);
 bone->AddElement(elC,0.278);
 bone->AddElement(elN,0.027);
 bone->AddElement(elO,0.410016);
 bone->AddElement(elMg,0.002);
 bone->AddElement(elP,0.07);
 bone->AddElement(elS,0.002);
 bone->AddElement(elCa,0.147);

 // Volumes

// NECK (our world Volume)

 G4double Rmin = 0.*cm;
 G4double Rmax = 6.*cm;
 G4double Dz = 8.*cm;
 G4double SPhi = 0.*deg;
 G4double DPhi = 360.*deg;


 G4Tubs* WaterNeck = new G4Tubs("Waterneck", Rmin, Rmax, Dz, SPhi, DPhi);
 G4LogicalVolume* WaterNeckLog = new
G4LogicalVolume(WaterNeck,matH2O,"WaterNeckLog",0,0,0); 
G4VPhysicalVolume* WaterNeckPhys = new
G4PVPlacement(0,G4ThreeVector(),"WaterNeckPhys",WaterNeckLog,0,false,0);
	
// TRACHEA (bone Tubs)

 G4double TRmin = 0.3*cm;
 G4double TRmax = 0.5*cm;
 G4double TDz = 7.*cm;
 G4double TSPhi = 0.*deg;
 G4double TDPhi = 360.*deg;

 G4Tubs* Trachea = new G4Tubs("Trachea",TRmin, TRmax, TDz, TSPhi, TDPhi);
 G4LogicalVolume* TracheaLog = new
 G4LogicalVolume(Trachea,bone,"TracheaLog",0,0,0); 
 G4VPhysicalVolume* TracheaPhys = new
 G4PVPlacement(0,G4ThreeVector(),TracheaLog,"TracheaPhys",WaterNeckLog,false,0);

 // Interno Trachea (Air Tubs)

  G4double ARmin = 0.*cm;
  G4double ARmax = 0.3*cm;
  G4double ADz = 7.*cm;
  G4double ASPhi = 0.*deg;
  G4double ADPhi = 360.*deg;
 //
  G4Tubs* IntTrachea = new G4Tubs("IntTrachea",ARmin, ARmax, ADz, ASPhi, ADPhi);
  G4LogicalVolume* IntTracheaLog = new
  G4LogicalVolume(IntTrachea,matAir,"IntTracheaLog",0,0,0); 
  G4VPhysicalVolume* IntTracheaPhys = new
  G4PVPlacement(0,G4ThreeVector(),IntTracheaLog,"IntTracheaPhys",TracheaLog,false,0);

 // Thyroid  (2 ellipticaltubes)
 
 // LeftThyroid (Tubo ellittico ruotato)
    G4double TlDx = 0.4*cm;
    G4double TlDy = 1.0*cm;
    G4double TlDz = 2.5*cm;
    

 G4EllipticalTube* LeftThyroid  = new G4EllipticalTube ("LeftThyroid", TlDx, TlDy, TlDz);
 G4RotationMatrix* rotD1 = new G4RotationMatrix();
  rotD1->rotateX (20.*deg);
 G4LogicalVolume*LeftThyroidLog = new
G4LogicalVolume(LeftThyroid,matSoftTissue,"LeftThyroidLog",0,0,0); 
 //G4VPhysicalVolume* LeftThyroidPhys = new
 //G4PVPlacement(rotD1,G4ThreeVector(0,-1.0*cm,0),LeftThyroidLog,"LeftThyroidPhys",WaterNeckLog,false,0);

 


 // RightThyroid (Tubo ellittico ruotato)
   
    G4double TrDx = 0.4*cm;
    G4double TrDy = 1.0*cm;
    G4double TrDz = 2.5*cm;

   G4RotationMatrix* rotD2 = new G4RotationMatrix();
    rotD2->rotateX (-20.*deg);

printf("\n tiroide\n");

 G4EllipticalTube* RightThyroid  = new G4EllipticalTube("RightThyroid", TrDx, TrDy, TrDz);
 G4LogicalVolume*RightThyroidLog = new
G4LogicalVolume(RightThyroid,matSoftTissue,"RightThyroidLog",0,0,0); 
 //G4VPhysicalVolume* RightThyroidPhys = new
 //G4PVPlacement(rotD2,G4ThreeVector(0,1.0*cm,0),RightThyroidLog,"RightThyroidPhys",WaterNeckLog,false,0);

// Now make an union of the two
 G4double deltay = TrDz*sin(40*deg);
 G4double deltaz = TrDz-(deltay/tan(40*deg));
 G4ThreeVector moveinunion( 0.0, deltay, deltaz );
 G4RotationMatrix* rotinunion = new G4RotationMatrix;
 rotinunion->rotateX(40*deg);
 G4UnionSolid* thyroidUnion = new G4UnionSolid( "unionThyroidSolid", LeftThyroid, RightThyroid,
                                                rotinunion, moveinunion ); 
 G4LogicalVolume* UnionThyroidLog = new G4LogicalVolume(thyroidUnion,matSoftTissue,"unionThyroidLog",0,0,0); 
 G4VPhysicalVolume* unionThyroidPhysical = new G4PVPlacement( rotD2, G4ThreeVector(0,0,0),UnionThyroidLog,"unionThyroidPhys",WaterNeckLog,false,0);
 
 G4ThreeVector moveinmother(0.7*cm,-0.58*deltay,0.0); 

 // Iodine Nodule
printf("\n nodulo1\n");

 G4Sphere* IodineNodule = new G4Sphere("IodineNodule",
0.*cm,0.5*cm,0.*deg,360.*deg,0.*deg,180.*deg);  G4LogicalVolume*
IodineNoduleLog = new
G4LogicalVolume(IodineNodule,matI131,"IodineNoduleLog",0,0,0); 
G4VPhysicalVolume* IodineNodulePhys = new
G4PVPlacement(0,G4ThreeVector(0,-1.0*cm,0),IodineNoduleLog,"IodineNodulePhys",LeftThyroidLog,false,0); 	  // Sensitive Detector and ReadOut geometry definition

printf("\n nodulo\n");

 G4SDManager* pSDManager = G4SDManager::GetSDMpointer();
 ThyroidSD* pThyroidSD = new ThyroidSD(m_SDName);
 if(pThyroidSD)
	{
	G4String ROGeometryName = "ThyroidROGeometry";
	ThyroidROGeometry* pThyroidROGeometry = new
ThyroidROGeometry(ROGeometryName,
m_theDx,m_theDy,
m_theDz); 
pThyroidROGeometry->BuildROGeometry(); 
pThyroidSD->SetROgeometry(pThyroidROGeometry); 
pSDManager->AddNewDetector(pThyroidSD);        
RightThyroidLog->SetSensitiveDetector(pThyroidSD); 
 }


 G4VisAttributes* IodineNoduleVisAtt= new G4VisAttributes(magenta);
   IodineNoduleVisAtt->SetVisibility(true);
   IodineNoduleVisAtt->SetForceWireframe(true);
  IodineNoduleLog->SetVisAttributes(IodineNoduleVisAtt);

 G4VisAttributes* LeftThyroidVisAtt= new G4VisAttributes(red);
 LeftThyroidVisAtt->SetVisibility(true);  
 LeftThyroidVisAtt->SetForceWireframe(true);
 LeftThyroidLog->SetVisAttributes( LeftThyroidVisAtt);

 G4VisAttributes* RightThyroidVisAtt= new G4VisAttributes(red);
 RightThyroidVisAtt->SetVisibility(true);  
 RightThyroidVisAtt->SetForceWireframe(true);
 RightThyroidLog->SetVisAttributes( RightThyroidVisAtt);

  G4VisAttributes* TracheaVisAtt= new G4VisAttributes(blue);
  TracheaVisAtt->SetVisibility(true);  
  TracheaVisAtt->SetForceWireframe(true);
  TracheaLog->SetVisAttributes( TracheaVisAtt);

  G4VisAttributes* IntTracheaVisAtt= new G4VisAttributes(grey);
  IntTracheaVisAtt->SetVisibility(true);  
  IntTracheaVisAtt->SetForceWireframe(true);
  IntTracheaLog->SetVisAttributes( IntTracheaVisAtt);

 return WaterNeckPhys;
}



















