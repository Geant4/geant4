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
// Rich advanced example for Geant4
// RichTbComponent.cc for Rich of LHCb
// History:
// Created: Sajan Easo (Sajan.Easo@cern.ch)
// Revision and changes: Patricia Mendez (Patricia.Mendez@cern.ch)
/////////////////////////////////////////////////////////////////////////////
#include <iostream>
#include "globals.hh"
#include "RichTbDetectorConstruction.hh"
#include "RichTbGeometryParameters.hh"
#include "RichTbMaterialParameters.hh"
#include "FilterTypeSpec.hh"
#include "RichTbComponent.hh"
#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4Sphere.hh"
#include "G4LogicalVolume.hh"
#include "G4VPhysicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4RotationMatrix.hh"
#include "G4ThreeVector.hh"
#include "G4Transform3D.hh"
#include "G4LogicalBorderSurface.hh"
#include "G4Material.hh"
#include "G4SubtractionSolid.hh"

RichTbComponent::RichTbComponent(){ ; }
RichTbComponent::RichTbComponent(RichTbMaterial* RMaterial, 
  RichTbHall* RTbHall , RichTbRunConfig* RConfig, 
  G4bool ConstructTrackingSwitch):
  RichTbAgelLVol(std::vector<G4LogicalVolume*>(MaxNumberOfAerogelTiles)),
  RichTbAgelPVol(std::vector<G4VPhysicalVolume*>(MaxNumberOfAerogelTiles)),
  RichTbAgelWrapTopLVol(std::vector<G4LogicalVolume*>(MaxNumberOfAerogelTiles)),
  RichTbAgelWrapTopPVol(std::vector<G4VPhysicalVolume*>(MaxNumberOfAerogelTiles)),
  RichTbAgelWrapBotLVol(std::vector<G4LogicalVolume*>(MaxNumberOfAerogelTiles)),
 RichTbAgelWrapBotPVol(std::vector<G4VPhysicalVolume*>(MaxNumberOfAerogelTiles))
 {  

    ConstructTrackingGeometrySwitch=ConstructTrackingSwitch;

  //First define the shapes.
  //Components of the Vessel.
  G4Tubs * VesselEnclosure
    = new G4Tubs("VesselEnclosure",VesselInnerRad,VesselOuterRad,
                 VesselHalfZ,VesselStartPhi,VesselDelPhi);
  G4RotationMatrix VesselBRot;  
  G4ThreeVector VesselBTrsl(VesselPosX,VesselPosY,VesselPosZ);

  G4Transform3D VesselBTransform(VesselBRot,VesselBTrsl);

  G4LogicalVolume* VesselEnclosureLog=
    new G4LogicalVolume(VesselEnclosure, RMaterial->getNitrogenGas(),
                        "VesselEnclosure",0,0,0);
 
  G4VPhysicalVolume* VesselEnclosurePhys=
    new G4PVPlacement(VesselBTransform,"VesselEnclosure",
		      VesselEnclosureLog,
                      RTbHall->getRichTbHallPhysicalVolume(),
                      false,0);
  //
  //
  if(ConstructTrackingSwitch){

  G4int NAgelTiles=RConfig-> GetNumberOfAerogelTiles();


  G4Box * RadFrame 
    = new G4Box("RadFrame",RadFrameHalfX,RadFrameHalfY,RadFrameHalfZ);

  G4Box * RadHoldUps 
    = new G4Box("RadHoldUps",RadHoldUpHalfX,RadHoldUpHalfY,RadHoldUpHalfZ);
  G4Tubs * RadHoldUpWin
    = new G4Tubs("RadHoldUpWin", RadWinUpInnerRad, RadWinUpOuterRad,
                  RadWinUpHalfZ,RadWinUpStartPhi,RadWinUpDelPhi);


  G4Box * RadHoldDns 
    = new G4Box("RadHoldDns",RadHoldDnHalfX,RadHoldDnHalfY,RadHoldDnHalfZ);
  G4Tubs * RadHoldDnWin
    = new G4Tubs("RadHoldDnWin", RadWinDnInnerRad, RadWinDnOuterRad,
                  RadWinDnHalfZ,RadWinDnStartPhi,RadWinDnDelPhi);


  G4int FilterNum = RConfig->GetFilterTNumber();
  G4Box * FilterBox=0;
  if(FilterNum >= 0 ) {
  G4double FilterHalfZ= FilterHalfZArray[FilterNum];
  FilterBox 
    = new G4Box("FilterBox", FilterHalfX,FilterHalfY,FilterHalfZ);

  }
  //Mirror
  G4double MirrThetaSize=2.0*std::asin(MirrorHorizontalChord/(MirrorRInner*2.0));
  G4double MirrThetaStart=halfpi*rad-MirrThetaSize/2.0;
  G4double MirrPhiSize=2.0*std::asin(MirrorVerticalChord/(MirrorRInner*2.0));
  G4double MirrPhiStart=-MirrPhiSize/2.0;

  G4Sphere* MirrorSphe= new G4Sphere("MirrorSphe",MirrorRInner,
				      MirrorROuter,MirrPhiStart,
                                      MirrPhiSize,MirrThetaStart,
                                      MirrThetaSize);

  //Now for the rotations and translations..

  G4RotationMatrix RadFrameRot;
  G4RotationMatrix RadWinUpRot,RadWinDnRot,RadHoldUpRot,RadHoldDnRot;

  G4ThreeVector RadFrameTrsl(RadFramePosX,RadFramePosY,RadFramePosZ);
  G4Transform3D RadFrameTransform(RadFrameRot,RadFrameTrsl);
  //
  //
  G4ThreeVector RadWinUpTrsl(RadWinUpShiftX,RadWinUpShiftY,RadWinUpShiftZ);
  G4Transform3D RadWinUpTransform(RadWinUpRot,RadWinUpTrsl);
  G4ThreeVector RadHoldUpTrsl(RadHoldUpPosX,RadHoldUpPosY,RadHoldUpPosZ);
  G4Transform3D RadHoldUpTransform(RadHoldUpRot,RadHoldUpTrsl);
  //
  G4ThreeVector RadWinDnTrsl(RadWinDnShiftX,RadWinDnShiftY,RadWinDnShiftZ);
  G4Transform3D RadWinDnTransform(RadWinDnRot,RadWinDnTrsl);
  G4ThreeVector RadHoldDnTrsl(RadHoldDnPosX,RadHoldDnPosY,RadHoldDnPosZ);
  G4Transform3D RadHoldDnTransform(RadHoldDnRot,RadHoldDnTrsl);
  //
  //
  //
  G4RotationMatrix MirrorRotationX, MirrorRotationY;
  G4double MirrorExtraTiltX = RConfig->getMirrorAddTiltX();
  G4double MirrorExtraTiltY = RConfig->getMirrorAddTiltY();
  G4double MirrorTotRotX= MirrorNominalRotX+MirrorExtraTiltX;
  G4double MirrorTotRotY= MirrorNominalRotY+MirrorExtraTiltY;
  G4double  MirrorPosX,  MirrorPosY,  MirrorPosZ;

  if( MirrorTotRotX !=0.0 ||  MirrorTotRotY != 0.0 ) {
     MirrorPosX =-VesselPosX-MirrorRInner*std::sin(MirrorTotRotY);
     MirrorPosY =-VesselPosY+MirrorRInner*std::sin(MirrorTotRotX);
     MirrorPosZ =MirrorShiftFromEnd-VesselHalfZ
      -MirrorRInner*(std::cos(MirrorTotRotY)+std::cos(MirrorTotRotX)-1.0);
  }else {
    MirrorPosX =-VesselPosX;
    MirrorPosY =-VesselPosY;
    MirrorPosZ =MirrorNominalPosZ;
  }
  G4ThreeVector MirrorPos(MirrorPosX,MirrorPosY,MirrorPosZ);
  if (MirrorTotRotX != 0.0 ) {
  MirrorRotationX.rotateX(MirrorTotRotX);
  }
  MirrorRotationY.rotateY(-halfpi*rad+MirrorTotRotY);
  G4Transform3D MirrorTransform(MirrorRotationX*MirrorRotationY,MirrorPos);

  //
       
  G4LogicalBorderSurface* VesselOuterSurface =
    new G4LogicalBorderSurface("VesselOuterSurface", 
			 RTbHall->getRichTbHallPhysicalVolume(),
                         VesselEnclosurePhys,
                         RMaterial->getOpticalEnclosureSurface());
  G4LogicalBorderSurface* VesselInnerSurface =
    new G4LogicalBorderSurface("VesselInnerSurface",
                         VesselEnclosurePhys, 
			 RTbHall->getRichTbHallPhysicalVolume(),
                         RMaterial->getOpticalEnclosureSurface());

  RichTbEnclosureOuterBSurf=VesselOuterSurface;
  RichTbEnclosureInnerBSurf=VesselInnerSurface;

  G4LogicalVolume* RadFrameLog=
    new G4LogicalVolume(RadFrame,RMaterial->getNitrogenGas(),
                        "RadFrame",0,0,0);

  G4VPhysicalVolume* RadFramePhys=
    new G4PVPlacement(RadFrameTransform,"RadFrame",RadFrameLog,
                     VesselEnclosurePhys,false,0);

  //Now for the holder upstream and downstream of the Aerogel Tiles.
  G4SubtractionSolid* RadUpW = 
    new G4SubtractionSolid("RadUpW",RadHoldUps,RadHoldUpWin,
			   RadWinUpTransform); 

  G4LogicalVolume* RadUpWLog=
    new G4LogicalVolume(RadUpW,RMaterial->getPlasticAg(),
                        "RadUpW",0,0,0);

  G4VPhysicalVolume* RadUpWPhys=
    new G4PVPlacement( RadHoldUpTransform,"RadUpW",RadUpWLog,
                     RadFramePhys ,false,0);
  //
  G4SubtractionSolid* RadDnW = 
    new G4SubtractionSolid("RadDnW",RadHoldDns,RadHoldDnWin,
			   RadWinDnTransform); 

  G4LogicalVolume* RadDnWLog=
    new G4LogicalVolume(RadDnW,RMaterial->getAluminium(),
                        "RadDnW",0,0,0);

  G4VPhysicalVolume* RadDnWPhys=
    new G4PVPlacement( RadHoldDnTransform,"RadDnW",RadDnWLog,
                     RadFramePhys ,false,0);
   

  //Now for the various aerogel tiles and the wraps above and
  //below them

  NumAerogelTiles= NAgelTiles; 
  G4RotationMatrix AgelRot,AgelWrapTopRot,AgelWrapBotRot;
  G4double AgelPosSZ[MaxNumberOfAerogelTiles];
  G4double AgelPosCurZ[MaxNumberOfAerogelTiles];
  for(G4int agNum=0; agNum<NAgelTiles; agNum++ ){
     G4int agTnumber= RConfig->GetCurAerogelTNumber(agNum);

     G4Box* Agel = new G4Box("Agel", AgelHalfX[agTnumber],
                           AgelHalfY[agTnumber],AgelHalfZ[agTnumber]);

     AgelPosSZ[agNum]=AgelHalfZ[agTnumber];

     AgelPosCurZ[agNum]=AgelEndZ-AgelHalfZ[agTnumber];
     for(G4int agt=0; agt<agNum; agt++ ) {

       AgelPosCurZ[agNum] -= (2*AgelPosSZ[agt]+ AgelTileGapZ);
     } 
     G4ThreeVector AgelTrsl(AgelPosX[agNum],
            AgelPosY[agNum],AgelPosCurZ[agNum]);
     G4Transform3D AgelTransform(AgelRot,AgelTrsl);

     G4Material* CurAgelMaterial= RMaterial->getAerogelMaterial(agTnumber);

     G4LogicalVolume* AgelLog=
       new G4LogicalVolume(Agel,CurAgelMaterial,"Agel",0,0,0);
       G4VPhysicalVolume* AgelPhys=
        new G4PVPlacement(AgelTransform,"Agel",AgelLog,
                     RadFramePhys,false,agNum);

     RichTbAgelLVol[agNum] =  AgelLog;
     RichTbAgelPVol[agNum] =  AgelPhys;

     G4Box* AgelWrapTop = new G4Box("AgelWrapTop", AgelWrapTopHalfX[agTnumber],
                           AgelWrapTopHalfY[agTnumber],
                           AgelWrapTopHalfZ[agTnumber]);
     G4Box* AgelWrapBot = new G4Box("AgelWrapBot", AgelWrapBotHalfX[agTnumber],
                       AgelWrapBotHalfY[agTnumber],
                       AgelWrapBotHalfZ[agTnumber]);

     G4double AgelWrapTopPosCurZ=AgelPosCurZ[agNum];
     G4double AgelWrapBotPosCurZ=AgelPosCurZ[agNum];
     G4double AgelWrapTopPosCurY=AgelHalfY[agTnumber]+ 
       AgelWrapTopHalfY[agTnumber];
     G4double AgelWrapBotPosCurY=-AgelHalfY[agTnumber]+ 
       -AgelWrapBotHalfY[agTnumber];
     G4ThreeVector AgelWrapTopTrsl(AgelWrapTopPosX[agNum],
                   AgelWrapTopPosCurY,AgelWrapTopPosCurZ);
     G4Transform3D AgelWrapTopTransform(AgelWrapTopRot,AgelWrapTopTrsl);
     G4ThreeVector AgelWrapBotTrsl(AgelWrapBotPosX[agNum],
                   AgelWrapBotPosCurY,AgelWrapBotPosCurZ);
     G4Transform3D AgelWrapBotTransform(AgelWrapBotRot,AgelWrapBotTrsl);

     G4LogicalVolume* AgelWrapTopLog=
        new G4LogicalVolume(AgelWrapTop, RMaterial->getPlasticAg(),
                        "AgelWrapTop",0,0,0);
      G4LogicalVolume* AgelWrapBotLog=
       new G4LogicalVolume(AgelWrapBot, RMaterial->getPlasticAg(),
                        "AgelWrapBot",0,0,0);

   G4VPhysicalVolume* AgelWrapTopPhys=
      new G4PVPlacement(AgelWrapTopTransform,"AgelWrapTop",
                    AgelWrapTopLog,RadFramePhys,false,agNum);

   G4VPhysicalVolume* AgelWrapBotPhys=
     new G4PVPlacement(AgelWrapBotTransform,"AgelWrapBot",
                  AgelWrapBotLog,RadFramePhys,false,agNum);

   RichTbAgelWrapTopLVol[agNum] =  AgelWrapTopLog;
   RichTbAgelWrapTopPVol[agNum] =  AgelWrapTopPhys;
   RichTbAgelWrapBotLVol[agNum] =  AgelWrapBotLog;
   RichTbAgelWrapBotPVol[agNum] =  AgelWrapBotPhys;

  }

  //  FilterType CurFil= RConfig->GetFilterType();
  G4int Filnum=  RConfig->GetFilterTNumber() ;
  G4LogicalVolume* FilterLog=0;
  G4VPhysicalVolume* FilterPhys=0;
  if(Filnum >= 0 ) {

  G4double FilterHalfZCur= FilterHalfZArray[Filnum];

  G4double FilterPosZ= FilterPosZNominal-FilterHalfZNominal+FilterHalfZCur;
  

  G4RotationMatrix FilterRot;

  G4ThreeVector FilterTrsl(FilterPosX,FilterPosY,FilterPosZ);
  G4Transform3D FilterTransform(FilterRot,FilterTrsl);

  G4Material*  CurFilterMaterial
    = RMaterial-> getRichTbFilterMaterial( Filnum);

     FilterLog=
      new G4LogicalVolume(FilterBox, CurFilterMaterial,
                        "FilterBox",0,0,0);

    FilterPhys=
    new G4PVPlacement(FilterTransform,"FilterBox",FilterLog,
                     RadFramePhys,false,0);

//  G4LogicalBorderSurface* FilterInnerSurface =
    new G4LogicalBorderSurface("RichTbFilterInnerSurface",
			       FilterPhys,VesselEnclosurePhys,
                               RMaterial->getOpticalFilterSurface());

//  G4LogicalBorderSurface* FilterOuterSurface =
    new G4LogicalBorderSurface("RichTbFilterOuterSurface",
			       VesselEnclosurePhys,FilterPhys,
                               RMaterial->getOpticalFilterSurface());

  }
  G4LogicalVolume* MirrorLog=
    new G4LogicalVolume(MirrorSphe,RMaterial->getMirrorQuartz(),
                        "MirrorSphe",0,0,0);


  G4VPhysicalVolume* MirrorPhys=
    new G4PVPlacement(MirrorTransform,"MirrorSphe",MirrorLog,
                     VesselEnclosurePhys,false,0);
  G4LogicalBorderSurface* MirrorSurface =
    new G4LogicalBorderSurface("RichTbMirrorSurface",
			       VesselEnclosurePhys,MirrorPhys,
                               RMaterial->getOpticalMirrorSurface());


  RichTbRadFrameLVol= RadFrameLog;
  RichTbRadFramePVol= RadFramePhys;

  RichTbRadUpWLVol= RadUpWLog;
  RichTbRadUpWPVol= RadUpWPhys;

  RichTbRadDnWLVol= RadDnWLog;
  RichTbRadDnWPVol= RadDnWPhys;

  if(Filnum >= 0 ) {
  RichTbFilterLVol =  FilterLog;
  RichTbFilterPVol =  FilterPhys;
  }
  RichTbMirrorLVol=MirrorLog;
  RichTbMirrorPVol=MirrorPhys;
  RichTbMirrorBSurf=MirrorSurface;
  }

  RichTbEnclosureLVol=VesselEnclosureLog;
  RichTbEnclosurePVol=VesselEnclosurePhys;


 }

RichTbComponent::~RichTbComponent() {; }









