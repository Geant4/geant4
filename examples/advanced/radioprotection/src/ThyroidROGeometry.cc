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
//    ************************************
//    *                                  *
//    *       ThyroidROGeometry.cc       *
//    *                                  *
//    ************************************

#include "ThyroidROGeometry.hh"

#include "G4LogicalVolume.hh"
#include "G4VPhysicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4PVReplica.hh"
#include "G4SDManager.hh"
#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4SubtractionSolid.hh"
#include "G4ThreeVector.hh"
#include "G4Material.hh"
#include "G4EllipticalTube.hh"
#include "G4RotationMatrix.hh"
#include "G4MaterialTable.hh"
#include "G4Element.hh"
#include "G4ElementTable.hh"
//....

ThyroidROGeometry::ThyroidROGeometry(G4String aString,G4double DetrDx, G4double DetrDy,G4double DetrDz)
  : G4VReadOutGeometry(aString),m_DetrDx(DetrDx),m_DetrDy(DetrDy),m_DetrDz(DetrDz)
{
}

//....

ThyroidROGeometry::~ThyroidROGeometry()
{
}

//....

G4VPhysicalVolume* ThyroidROGeometry::Build()
{
 // A dummy material is used to fill the volumes of the readout geometry.
 // (It will be allowed to set a NULL pointer in volumes of such virtual
 // division in future, since this material is irrelevant for tracking.)

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
   
 //SoftTissue



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



 
// World (Neck)

 G4double Rmin = 0.*cm;
 G4double Rmax = 6.*cm;
 G4double Dz = 8.*cm;
 G4double SPhi = 0.*deg;
 G4double DPhi = 360.*deg;

 G4Tubs *ROWaterNeck = new G4Tubs("ROWaterNeck", Rmin, Rmax, Dz, SPhi, DPhi);
 G4LogicalVolume *ROWaterNeckLog = new
G4LogicalVolume(ROWaterNeck,matH2O,"ROWaterNeckLog",0,0,0);  G4VPhysicalVolume
*ROWaterNeckPhys = new
G4PVPlacement(0,G4ThreeVector(),"ROWaterNeckPhys",ROWaterNeckLog,0,false,0);

 // Sensible Thyroid 


 G4EllipticalTube *RODetector = new G4EllipticalTube("RODetector", m_DetrDx,m_DetrDy, m_DetrDz);
 G4RotationMatrix* rotD2 = new G4RotationMatrix ();
 rotD2->rotateX ( -20.*deg);
 G4LogicalVolume *RODetectorLog = new
G4LogicalVolume(RODetector,matSoftTissue,"RODetectorLog",0,0,0); 
G4VPhysicalVolume *RODetectorPhys = new
G4PVPlacement(rotD2,G4ThreeVector(0,1,0),"DetectorPhys",RODetectorLog,ROWaterNeckPhys,false,0);

  


 return ROWaterNeckPhys;
}










