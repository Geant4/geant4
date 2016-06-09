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
// $Id: MedLinacMLCDecorator.cc,v 1.3 2006/06/29 16:04:27 gunter Exp $
//
// Code developed by: M. Piergentili
//
//
#include "MedLinacVGeometryComponent.hh"
#include "G4Material.hh"
#include "MedLinacMLCDecorator.hh"
#include "MedLinacDecorator.hh"
#include "MedLinacMLCMessenger.hh"

#include "globals.hh"
#include "G4LogicalVolume.hh"
#include "G4RotationMatrix.hh"
#include "G4Transform3D.hh"
#include "G4VPhysicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4Material.hh"
#include "G4MaterialPropertiesTable.hh"
#include "G4MaterialTable.hh"
#include "G4MaterialPropertyVector.hh"
#include "G4Element.hh"
#include "G4ElementTable.hh"
#include "G4Box.hh"
#include "G4Cons.hh"
#include "G4Tubs.hh"
#include "G4ThreeVector.hh"
#include "G4VisAttributes.hh"
#include "G4GeometryManager.hh"
#include "G4BooleanSolid.hh"
#include "G4SubtractionSolid.hh"
#include "G4UnionSolid.hh"
#include "G4VSolid.hh"
#include "G4PhysicalVolumeStore.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4SolidStore.hh"
#include "G4Colour.hh"

#include "G4ios.hh"
MedLinacMLCDecorator::MedLinacMLCDecorator(MedLinacVGeometryComponent* comp)
  : MedLinacDecorator(comp),
    leafLog(0), 
    leafAPhys(0), leafBPhys(0)
{

   // default parameter values

  a1y = 40.*cm;
  a2y = 40.*cm;
  a3y = 40.*cm;
  a4y = 40.*cm;
  a5y = 40.*cm;
  a6y = 40.*cm;
  a7y = 40.*cm;
  a8y = 40.*cm;
  a9y = 40.*cm;
  a10y = 40.*cm;
  a11y = 40.*cm;
  a12y = 40.*cm;
  a13y = 40.*cm;
  a14y = 40.*cm;
  a15y = 40.*cm;
  a16y = 40.*cm;
  a17y = 40.*cm;
  a18y = 40.*cm;
  a19y = 40.*cm;
  a20y = 40.*cm;
  a21y = 40.*cm;
  a22y = 40.*cm;
  a23y = 40.*cm;
  a24y = 40.*cm;
  a25y = 40.*cm;
  a26y = 40.*cm;
  a27y = 40.*cm;
  a28y = 40.*cm;
  a29y = 40.*cm;
  a30y = 40.*cm;
  a31y = 40.*cm;
  a32y = 40.*cm;
  a33y = 40.*cm;
  a34y = 40.*cm;
  a35y = 40.*cm;
  a36y = 40.*cm;
  a37y = 40.*cm;
  a38y = 40.*cm;
  a39y = 40.*cm;
  a40y = 40.*cm;

  b1y = 40.*cm;
  b2y = 40.*cm;
  b3y = 40.*cm;
  b4y = 40.*cm;
  b5y = 40.*cm;
  b6y = 40.*cm;
  b7y = 40.*cm;
  b8y = 40.*cm;
  b9y = 40.*cm;
  b10y = 40.*cm;
  b11y = 40.*cm;
  b12y = 40.*cm;
  b13y = 40.*cm;
  b14y = 40.*cm;
  b15y = 40.*cm;
  b16y = 40.*cm;
  b17y = 40.*cm;
  b18y = 40.*cm;
  b19y = 40.*cm;
  b20y = 40.*cm;
  b21y = 40.*cm;
  b22y = 40.*cm;
  b23y = 40.*cm;
  b24y = 40.*cm;
  b25y = 40.*cm;
  b26y = 40.*cm;
  b27y = 40.*cm;
  b28y = 40.*cm;
  b29y = 40.*cm;
  b30y = 40.*cm;
  b31y = 40.*cm;
  b32y = 40.*cm;
  b33y = 40.*cm;
  b34y = 40.*cm;
  b35y = 40.*cm;
  b36y = 40.*cm;
  b37y = 40.*cm;
  b38y = 40.*cm;
  b39y = 40.*cm;
  b40y = 40.*cm;

  MLCMessenger = new MedLinacMLCMessenger(this);
}
MedLinacMLCDecorator::~MedLinacMLCDecorator()
{
  delete MLCMessenger;
    ;
}
void MedLinacMLCDecorator::ConstructComponent(G4VPhysicalVolume* world, G4VPhysicalVolume* vacuumBlock)
{
   MedLinacDecorator::ConstructComponent(world,vacuumBlock);
   ConstructMultiLeafCollimator(world,vacuumBlock);
}

void MedLinacMLCDecorator::DestroyComponent()
{
  ;
}
void MedLinacMLCDecorator::ConstructMultiLeafCollimator(G4VPhysicalVolume* world, G4VPhysicalVolume*)
{

  //    materials

  G4double a;  // atomic mass
  G4double z;  // atomic number
  G4double density;
  G4String name;
 
   density = 18.*g/cm3;
  a = 183.85*g/mole;
  G4Material* W = new G4Material(name="Tungsten" , z=74., a, density);

  //    colors
  G4Colour  cyan    (0.0, 1.0, 1.0);

 //---------rotation matrix leaf end--------

  G4RotationMatrix*  rotateLeaf=new G4RotationMatrix();
  rotateLeaf->rotateY(90.0*deg);

  //---------rotation matrix leaves B--------

  G4RotationMatrix*  rotateLeavesB=new G4RotationMatrix();
  rotateLeavesB->rotateX(180.0*deg);
 
   //    volumes
  //    beam line along z axis
//----------------------

  //     single leaf

  G4double preLeafDim_x = 2.715*mm;
  G4double preLeafDim_y = 64.72066*mm;
  G4double preLeafDim_z = 29.3*mm;
  G4Box* preLeaf_box = new G4Box("preLeaf_box",preLeafDim_x,preLeafDim_y,preLeafDim_z);
  
  G4double preLeafPos_x = 0.0*m;
  G4double preLeafPos_y = 0.0*m;
  G4double preLeafPos_z = -50.*cm;
  G4DisplacedSolid* disPreLeaf = new G4DisplacedSolid("disPreLeaf",preLeaf_box,0,
                                 G4ThreeVector(preLeafPos_x,preLeafPos_y,preLeafPos_z));
  
  G4double innerRadiusOfTheLeafEnd = 0.01*mm;
  G4double outerRadiusOfTheLeafEnd = 80.*mm;
  G4double hightOfTheLeafEnd = 2.715*mm;
  G4double startAngleOfTheLeafEnd = 0.*deg;
  G4double spanningAngleOfTheLeafEnd = 180.*deg;
  G4Tubs* aLeafEnd = new G4Tubs("aLeafEnd",innerRadiusOfTheLeafEnd,
                                    outerRadiusOfTheLeafEnd,hightOfTheLeafEnd,
				    startAngleOfTheLeafEnd,spanningAngleOfTheLeafEnd);
    
  G4double gapDim_x = 90.00*mm;
  G4double gapDim_y = 74.18*mm;
  G4double gapDim_z = 3.0*mm;
  G4Box* gap = new G4Box("gap",gapDim_x,gapDim_y,gapDim_z);
  
  G4SubtractionSolid* leafEnd = new G4SubtractionSolid("leafEnd",aLeafEnd,gap);
  
  G4double leafEndPosX = 0.*cm;
  G4double leafEndPosY = -9.4413249*mm;
  G4double leafEndPosZ = -50.0*cm;
    
  G4DisplacedSolid* disLeafEnd = new G4DisplacedSolid("disLeafEnd",leafEnd,rotateLeaf,
                                 G4ThreeVector(leafEndPosX,leafEndPosY,leafEndPosZ)) ;

  G4UnionSolid* fullLeaf = new G4UnionSolid("fullLeaf", disPreLeaf, disLeafEnd);


  //++++++++++++++++++++++++++++

  G4double cutADim_x = 0.20*mm;
  G4double cutADim_y = 68.0*mm;
  G4double cutADim_z = 15.925*mm;
  G4Box* cutA = new G4Box("cutA",cutADim_x,cutADim_y,cutADim_z);

  G4double cutAPosX = -2.515*mm;
  G4double cutAPosY = 3.0*mm;
  G4double cutAPosZ = -48.6625*cm;
  G4DisplacedSolid* disCutA = new G4DisplacedSolid("disCutA",cutA,0,
                                 G4ThreeVector(cutAPosX,cutAPosY,cutAPosZ)) ;

  G4SubtractionSolid* nearlyLeaf = new G4SubtractionSolid("nearlyLeaf",fullLeaf,disCutA);
  //++++++++++++++++++++++++++++
  G4double cutBDim_x = 0.1625*mm;
  G4double cutBDim_y = 68.0*mm;
  G4double cutBDim_z = 14.415*mm;
  G4Box* cutB = new G4Box("cutB",cutBDim_x,cutBDim_y,cutBDim_z);

  G4double cutBPosX = 2.5525*mm;
  G4double cutBPosY = 3.0*mm;
  G4double cutBPosZ = -51.4885*cm;
  G4DisplacedSolid* disCutB = new G4DisplacedSolid("disCutB",cutB,0,
                                 G4ThreeVector(cutBPosX,cutBPosY,cutBPosZ)) ;

  G4SubtractionSolid* leaf = new G4SubtractionSolid("leaf",nearlyLeaf,disCutB);

  //++++++++++++++++++++++++++++



  leafLog = new G4LogicalVolume(leaf,W,"leafLog",0,0,0);

  G4double leafAPosX[41];
  G4double leafAPosY[41];
  G4double leafAPosYF[41];

  leafAPosX[0]= -111.335 *mm;

  leafAPosY[0]= 40.*cm;
  leafAPosY[1]= a1y;
  leafAPosY[2]= a2y;
  leafAPosY[3]= a3y;
  leafAPosY[4]= a4y;
  leafAPosY[5]= a5y;
  leafAPosY[6]= a6y;
  leafAPosY[7]= a7y;
  leafAPosY[8]= a8y;
  leafAPosY[9]= a9y;
  leafAPosY[10]= a10y;
  leafAPosY[11]= a11y;
  leafAPosY[12]= a12y;
  leafAPosY[13]= a13y;
  leafAPosY[14]= a14y;
  leafAPosY[15]= a15y;
  leafAPosY[16]= a16y;
  leafAPosY[17]= a17y;
  leafAPosY[18]= a18y;
  leafAPosY[19]= a19y;
  leafAPosY[20]= a20y;
  leafAPosY[21]= a21y;
  leafAPosY[22]= a22y;
  leafAPosY[23]= a23y;
  leafAPosY[24]= a24y;
  leafAPosY[25]= a25y;
  leafAPosY[26]= a26y;
  leafAPosY[27]= a27y;
  leafAPosY[28]= a28y;
  leafAPosY[29]= a29y;
  leafAPosY[30]= a30y;
  leafAPosY[31]= a31y;
  leafAPosY[32]= a32y;
  leafAPosY[33]= a33y;
  leafAPosY[34]= a34y;
  leafAPosY[35]= a35y;
  leafAPosY[36]= a36y;
  leafAPosY[37]= a37y;
  leafAPosY[38]= a38y;
  leafAPosY[39]= a39y;
  leafAPosY[40]= a40y;


  G4int CopyNbA = 0;

  for ( G4int i=1; i <= 40 ; i++ )
  { 
  leafAPosYF[i] = -(((leafAPosY[i]*(100.-47.25))/100.)+70.5587*mm);
  
  G4double leafAPosZ = 97.25 *cm;
    leafAPosX[i]=leafAPosX[i-1]+5.431 *mm;
    leafAPhys = new G4PVPlacement(0,G4ThreeVector(leafAPosX[i],leafAPosYF[i],leafAPosZ),
  		       		       "leafA",leafLog,world,false,CopyNbA);
    CopyNbA = CopyNbA+1; 
  }

  G4double leafBPosX[41];
  G4double leafBPosY[41];
  G4double leafBPosYF[41];

  G4double leafBPosZ = -2.75 *cm;
  leafBPosX[0]= -111.335*mm;
  G4int CopyNbB = 0;


  leafBPosY[0]= 40.*cm;
  leafBPosY[1]= b1y;
  leafBPosY[2]= b2y;
  leafBPosY[3]= b3y;
  leafBPosY[4]= b4y;
  leafBPosY[5]= b5y;
  leafBPosY[6]= b6y;
  leafBPosY[7]= b7y;
  leafBPosY[8]= b8y;
  leafBPosY[9]= b9y;
  leafBPosY[10]= b10y;
  leafBPosY[11]= b11y;
  leafBPosY[12]= b12y;
  leafBPosY[13]= b13y;
  leafBPosY[14]= b14y;
  leafBPosY[15]= b15y;
  leafBPosY[16]= b16y;
  leafBPosY[17]= b17y;
  leafBPosY[18]= b18y;
  leafBPosY[19]= b19y;
  leafBPosY[20]= b20y;
  leafBPosY[21]= b21y;
  leafBPosY[22]= b22y;
  leafBPosY[23]= b23y;
  leafBPosY[24]= b24y;
  leafBPosY[25]= b25y;
  leafBPosY[26]= b26y;
  leafBPosY[27]= b27y;
  leafBPosY[28]= b28y;
  leafBPosY[29]= b29y;
  leafBPosY[30]= b30y;
  leafBPosY[31]= b31y;
  leafBPosY[32]= b32y;
  leafBPosY[33]= b33y;
  leafBPosY[34]= b34y;
  leafBPosY[35]= b35y;
  leafBPosY[36]= b36y;
  leafBPosY[37]= b37y;
  leafBPosY[38]= b38y;
  leafBPosY[39]= b39y;
  leafBPosY[40]= b40y;
 

  for ( G4int i=1; i <= 40 ; i++ )
  {
  leafBPosYF[i] = (((leafBPosY[i]*(100.-47.25))/100.)+70.5587*mm);
  leafBPosX[i]=leafBPosX[i-1]+5.431 *mm;
  leafBPhys = new G4PVPlacement(rotateLeavesB,
  			    G4ThreeVector(leafBPosX[i],leafBPosYF[i],leafBPosZ),
  		       		       "leafB",leafLog,world,false,CopyNbB);
    CopyNbB = CopyNbB+1; 
  }


  PrintParametersMLC(); 

  //    Visualization attributes 

  G4VisAttributes* simpleTungstenSVisAtt= new G4VisAttributes(cyan);
  simpleTungstenSVisAtt->SetVisibility(true);
  simpleTungstenSVisAtt->SetForceSolid(true);
  leafLog->SetVisAttributes(simpleTungstenSVisAtt);

 }
void MedLinacMLCDecorator::PrintParametersMLC()
{
  G4cout <<"leaf A1 position "<< a1y/cm << " cm "<<G4endl ;  
  G4cout <<"leaf A2 position "<< a2y/cm << " cm "<<G4endl ; 
  G4cout <<"leaf A3 position "<< a3y/cm << " cm "<<G4endl ;  
  G4cout <<"leaf A4 position "<< a4y/cm << " cm "<<G4endl ; 
  G4cout <<"leaf A5 position "<< a5y/cm << " cm "<<G4endl ;  
  G4cout <<"leaf A6 position "<< a6y/cm << " cm "<<G4endl ; 
  G4cout <<"leaf A7 position "<< a7y/cm << " cm "<<G4endl ;  
  G4cout <<"leaf A8 position "<< a8y/cm << " cm "<<G4endl ; 
  G4cout <<"leaf A9 position "<< a9y/cm << " cm "<<G4endl ;  
  G4cout <<"leaf A10 position "<< a10y/cm << " cm "<<G4endl ; 
  G4cout <<"leaf A11 position "<< a11y/cm << " cm "<<G4endl ;  
  G4cout <<"leaf A12 position "<< a12y/cm << " cm "<<G4endl ; 
  G4cout <<"leaf A13 position "<< a13y/cm << " cm "<<G4endl ;  
  G4cout <<"leaf A14 position "<< a14y/cm << " cm "<<G4endl ; 
  G4cout <<"leaf A15 position "<< a15y/cm << " cm "<<G4endl ;  
  G4cout <<"leaf A16 position "<< a16y/cm << " cm "<<G4endl ; 
  G4cout <<"leaf A17 position "<< a17y/cm << " cm "<<G4endl ;  
  G4cout <<"leaf A18 position "<< a18y/cm << " cm "<<G4endl ; 
  G4cout <<"leaf A19 position "<< a19y/cm << " cm "<<G4endl ;  
  G4cout <<"leaf A20 position "<< a20y/cm << " cm "<<G4endl ; 
  G4cout <<"leaf A21 position "<< a21y/cm << " cm "<<G4endl ;  
  G4cout <<"leaf A22 position "<< a22y/cm << " cm "<<G4endl ; 
  G4cout <<"leaf A23 position "<< a23y/cm << " cm "<<G4endl ;  
  G4cout <<"leaf A24 position "<< a24y/cm << " cm "<<G4endl ; 
  G4cout <<"leaf A25 position "<< a25y/cm << " cm "<<G4endl ;  
  G4cout <<"leaf A26 position "<< a26y/cm << " cm "<<G4endl ; 
  G4cout <<"leaf A27 position "<< a27y/cm << " cm "<<G4endl ;  
  G4cout <<"leaf A28 position "<< a28y/cm << " cm "<<G4endl ; 
  G4cout <<"leaf A29 position "<< a29y/cm << " cm "<<G4endl ;  
  G4cout <<"leaf A30 position "<< a30y/cm << " cm "<<G4endl ;
  G4cout <<"leaf A31 position "<< a31y/cm << " cm "<<G4endl ;  
  G4cout <<"leaf A32 position "<< a32y/cm << " cm "<<G4endl ; 
  G4cout <<"leaf A33 position "<< a33y/cm << " cm "<<G4endl ;  
  G4cout <<"leaf A34 position "<< a34y/cm << " cm "<<G4endl ; 
  G4cout <<"leaf A35 position "<< a35y/cm << " cm "<<G4endl ;  
  G4cout <<"leaf A36 position "<< a36y/cm << " cm "<<G4endl ; 
  G4cout <<"leaf A37 position "<< a37y/cm << " cm "<<G4endl ;  
  G4cout <<"leaf A38 position "<< a38y/cm << " cm "<<G4endl ; 
  G4cout <<"leaf A39 position "<< a39y/cm << " cm "<<G4endl ;  
  G4cout <<"leaf A40 position "<< a40y/cm << " cm "<<G4endl ;

  G4cout <<"leaf B1 position "<< b1y/cm << " cm "<<G4endl ;  
  G4cout <<"leaf B2 position "<< b2y/cm << " cm "<<G4endl ; 
  G4cout <<"leaf B3 position "<< b3y/cm << " cm "<<G4endl ;  
  G4cout <<"leaf B4 position "<< b4y/cm << " cm "<<G4endl ; 
  G4cout <<"leaf B5 position "<< b5y/cm << " cm "<<G4endl ;  
  G4cout <<"leaf B6 position "<< b6y/cm << " cm "<<G4endl ; 
  G4cout <<"leaf B7 position "<< b7y/cm << " cm "<<G4endl ;  
  G4cout <<"leaf B8 position "<< b8y/cm << " cm "<<G4endl ; 
  G4cout <<"leaf B9 position "<< b9y/cm << " cm "<<G4endl ;  
  G4cout <<"leaf B10 position "<< b10y/cm << " cm "<<G4endl ; 
  G4cout <<"leaf B11 position "<< b11y/cm << " cm "<<G4endl ;  
  G4cout <<"leaf B12 position "<< b12y/cm << " cm "<<G4endl ; 
  G4cout <<"leaf B13 position "<< b13y/cm << " cm "<<G4endl ;  
  G4cout <<"leaf B14 position "<< b14y/cm << " cm "<<G4endl ; 
  G4cout <<"leaf B15 position "<< b15y/cm << " cm "<<G4endl ;  
  G4cout <<"leaf B16 position "<< b16y/cm << " cm "<<G4endl ; 
  G4cout <<"leaf B17 position "<< b17y/cm << " cm "<<G4endl ;  
  G4cout <<"leaf B18 position "<< b18y/cm << " cm "<<G4endl ; 
  G4cout <<"leaf B19 position "<< b19y/cm << " cm "<<G4endl ;  
  G4cout <<"leaf B20 position "<< b20y/cm << " cm "<<G4endl ; 
  G4cout <<"leaf B21 position "<< b21y/cm << " cm "<<G4endl ;  
  G4cout <<"leaf B22 position "<< b22y/cm << " cm "<<G4endl ; 
  G4cout <<"leaf B23 position "<< b23y/cm << " cm "<<G4endl ;  
  G4cout <<"leaf B24 position "<< b24y/cm << " cm "<<G4endl ; 
  G4cout <<"leaf B25 position "<< b25y/cm << " cm "<<G4endl ;  
  G4cout <<"leaf B26 position "<< b26y/cm << " cm "<<G4endl ; 
  G4cout <<"leaf B27 position "<< b27y/cm << " cm "<<G4endl ;  
  G4cout <<"leaf B28 position "<< b28y/cm << " cm "<<G4endl ; 
  G4cout <<"leaf B29 position "<< b29y/cm << " cm "<<G4endl ;  
  G4cout <<"leaf B30 position "<< b30y/cm << " cm "<<G4endl ;
  G4cout <<"leaf B31 position "<< b31y/cm << " cm "<<G4endl ;  
  G4cout <<"leaf B32 position "<< b32y/cm << " cm "<<G4endl ; 
  G4cout <<"leaf B33 position "<< b33y/cm << " cm "<<G4endl ;  
  G4cout <<"leaf B34 position "<< b34y/cm << " cm "<<G4endl ; 
  G4cout <<"leaf B35 position "<< b35y/cm << " cm "<<G4endl ;  
  G4cout <<"leaf B36 position "<< b36y/cm << " cm "<<G4endl ; 
  G4cout <<"leaf B37 position "<< b37y/cm << " cm "<<G4endl ;  
  G4cout <<"leaf B38 position "<< b38y/cm << " cm "<<G4endl ; 
  G4cout <<"leaf B39 position "<< b39y/cm << " cm "<<G4endl ;  
  G4cout <<"leaf B40 position "<< b40y/cm << " cm "<<G4endl ;
}
void MedLinacMLCDecorator::SetLeafName (G4String newleaf_name)
{ 
  leaf_name= newleaf_name;
}
void MedLinacMLCDecorator::SetPos_y (G4double pos)
{ 
  if (leaf_name=="a1")
    a1y = pos;
  if (leaf_name=="a2")
    a2y = pos;
  if (leaf_name=="a3")
    a3y = pos;
  if (leaf_name=="a4")
    a4y = pos;
  if (leaf_name=="a5")
    a5y = pos;
  if (leaf_name=="a6")
    a6y = pos;
  if (leaf_name=="a7")
    a7y = pos;
  if (leaf_name=="a8")
    a8y = pos;
  if (leaf_name=="a9")
    a9y = pos;
  if (leaf_name=="a10")
    a10y = pos;
  if (leaf_name=="a11")
    a11y = pos;
  if (leaf_name=="a12")
    a12y = pos;
  if (leaf_name=="a13")
    a13y = pos;
  if (leaf_name=="a14")
    a14y = pos;
  if (leaf_name=="a15")
    a15y = pos;
  if (leaf_name=="a16")
    a16y = pos;
  if (leaf_name=="a17")
    a17y = pos;
  if (leaf_name=="a18")
    a18y = pos;
  if (leaf_name=="a19")
    a19y = pos;
  if (leaf_name=="a20")
    a20y = pos;
  if (leaf_name=="a21")
    a21y = pos;
  if (leaf_name=="a22")
    a22y = pos;
  if (leaf_name=="a23")
    a23y = pos;
  if (leaf_name=="a24")
    a24y = pos;
  if (leaf_name=="a25")
    a25y = pos;
  if (leaf_name=="a26")
    a26y = pos;
  if (leaf_name=="a27")
    a27y = pos;
  if (leaf_name=="a28")
    a28y = pos;
  if (leaf_name=="a29")
    a29y = pos;
  if (leaf_name=="a30")
    a30y = pos;
  if (leaf_name=="a31")
    a31y = pos;
  if (leaf_name=="a32")
    a32y = pos;
  if (leaf_name=="a33")
    a33y = pos;
  if (leaf_name=="a34")
    a34y = pos;
  if (leaf_name=="a35")
    a35y = pos;
  if (leaf_name=="a36")
    a36y = pos;
  if (leaf_name=="a37")
    a37y = pos;
  if (leaf_name=="a38")
    a38y = pos;
  if (leaf_name=="a39")
    a39y = pos;
  if (leaf_name=="a40")
    a40y = pos;




  if (leaf_name=="b1")
    b1y = pos;
  if (leaf_name=="b2")
    b2y = pos;
  if (leaf_name=="b3")
    b3y = pos;
  if (leaf_name=="b4")
    b4y = pos;
  if (leaf_name=="b5")
    b5y = pos;
  if (leaf_name=="b6")
    b6y = pos;
  if (leaf_name=="b7")
    b7y = pos;
  if (leaf_name=="b8")
    b8y = pos;
  if (leaf_name=="b9")
    b9y = pos;
  if (leaf_name=="b10")
    b10y = pos;
  if (leaf_name=="b11")
    b11y = pos;
  if (leaf_name=="b12")
    b12y = pos;
  if (leaf_name=="b13")
    b13y = pos;
  if (leaf_name=="b14")
    b14y = pos;
  if (leaf_name=="b15")
    b15y = pos;
  if (leaf_name=="b16")
    b16y = pos;
  if (leaf_name=="b17")
    b17y = pos;
  if (leaf_name=="b18")
    b18y = pos;
  if (leaf_name=="b19")
    b19y = pos;
  if (leaf_name=="b20")
    b20y = pos;
  if (leaf_name=="b21")
    b21y = pos;
  if (leaf_name=="b22")
    b22y = pos;
  if (leaf_name=="b23")
    b23y = pos;
  if (leaf_name=="b24")
    b24y = pos;
  if (leaf_name=="b25")
    b25y = pos;
  if (leaf_name=="b26")
    b26y = pos;
  if (leaf_name=="b27")
    b27y = pos;
  if (leaf_name=="b28")
    b28y = pos;
  if (leaf_name=="b29")
    b29y = pos;
  if (leaf_name=="b30")
    b30y = pos;
  if (leaf_name=="b31")
    b31y = pos;
  if (leaf_name=="b32")
    b32y = pos;
  if (leaf_name=="b33")
    b33y = pos;
  if (leaf_name=="b34")
    b34y = pos;
  if (leaf_name=="b35")
    b35y = pos;
  if (leaf_name=="b36")
    b36y = pos;
  if (leaf_name=="b37")
    b37y = pos;
  if (leaf_name=="b38")
    b38y = pos;
  if (leaf_name=="b39")
    b39y = pos;
  if (leaf_name=="b40")
    b40y = pos;
}
