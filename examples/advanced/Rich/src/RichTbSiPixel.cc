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
// Rich advanced example for Geant4
// RichTbSiPixel.cc for Rich of LHCb
// History:
// Created: Sajan Easo (Sajan.Easo@cern.ch)
// Revision and changes: Patricia Mendez (Patricia.Mendez@cern.ch)
/////////////////////////////////////////////////////////////////////////////
#include <iostream>
#include "globals.hh"
#include "G4Box.hh"
#include "G4Trd.hh"
#include "G4RotationMatrix.hh"
#include "G4ThreeVector.hh"
#include "G4Transform3D.hh"
#include "G4LogicalVolume.hh"
#include "G4VPhysicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4UnionSolid.hh"
#include "RichTbSiPixel.hh"
#include "G4SDManager.hh"
#include "RichTbSD.hh"
#include "RichTbRODummySD.hh"
#include "RichTbMaterialParameters.hh"
#include "G4IntersectionSolid.hh"
RichTbSiPixel::RichTbSiPixel() {;}
RichTbSiPixel::~RichTbSiPixel() {;}
RichTbSiPixel::RichTbSiPixel(RichTbMaterial* RMaterial,
                         G4VPhysicalVolume* MotherOfSiPixel,
                         G4int IHpdNum, G4int Isector ,G4int ipixel){


  //
  //
  ICurSectorNumber=Isector;
  HpdSiPixelnum=ipixel;
  ICurHpdNumber=IHpdNum;


    G4Trd* HpdSiSectTrapInt=
      new G4Trd("HpdSiSectTrapInt",SiSectTrapHalfX1, SiSectTrapHalfX2,
		SiSectTrapHalfY1,SiSectTrapHalfY2,SiSectTrapHalfZ);

    G4Box* HpdSiPxBox;
    if(ipixel ==  BigPixelNum ) {
    HpdSiPxBox=new G4Box("HpdSiPxBigBox",
                         SiPixelHalfX,SiPixelHalfY,SiBigPixelHalfZ);
    }else {
    HpdSiPxBox=new G4Box("HpdSiPxBox",SiPixelHalfX,SiPixelHalfY,SiPixelHalfZ);
    }
  //
  //For the rotation and translation of the Si Pixel.
    G4RotationMatrix HpdSiPxRotation, HpdSiSectTrapIntRotation;  
    G4double RichTbHpdSiPxShiftY=HpdSiPixPosZ;
    G4double RichTbHpdSiPxShiftX
          =RichTbPadHpdSiPixPos(ipixel).getPadHpdSiPixPosX();
    G4double RichTbHpdSiPxShiftZ
          =RichTbPadHpdSiPixPos(ipixel).getPadHpdSiPixPosY()-SiSectTrapHalfZ;

    G4ThreeVector HpdSiPxTrsl( RichTbHpdSiPxShiftX,
                               RichTbHpdSiPxShiftY,
                               RichTbHpdSiPxShiftZ);

    G4Transform3D HpdSiPxTransform(HpdSiPxRotation,HpdSiPxTrsl);

    G4ThreeVector HpdSiSectTrapIntTrsl(0.0*mm, 0.0*mm, 0.0*mm);
    G4Transform3D HpdSiSectTrapIntTransform( HpdSiSectTrapIntRotation,
					     HpdSiSectTrapIntTrsl);
    //Now for the placements

    G4LogicalVolume*  HpdSiPxLog;
    G4VPhysicalVolume* HpdSiPxPhys;
    // The pixels at the pixel edge are intersected with the
    // sector. The BigPixel is also at the edge.

    if(PixelAtSectEdge[ipixel]){
    G4IntersectionSolid* HpdSiPxTrapBox
         = new G4IntersectionSolid("PixelSectTrapInt", HpdSiSectTrapInt,
			HpdSiPxBox, HpdSiPxTransform);
    

    HpdSiPxLog=
      new  G4LogicalVolume( HpdSiPxTrapBox,RMaterial->getHpdSiDetMaterial(),
			  "PixelSectTrapInt",0,0,0);         

    HpdSiPxPhys=
           new G4PVPlacement(HpdSiSectTrapIntTransform,"PixelSectTrapInt", 
                HpdSiPxLog,MotherOfSiPixel,false,ipixel);

    }else{

    HpdSiPxLog=
      new  G4LogicalVolume( HpdSiPxBox,RMaterial->getHpdSiDetMaterial(),
			  "HpdSiPxBox",0,0,0);         

    HpdSiPxPhys=
           new G4PVPlacement(HpdSiPxTransform,"HpdSiPxBox", 
                HpdSiPxLog,MotherOfSiPixel,false,ipixel);

    }
    //Now for the readout dummy sensitive detector
    // This just flags the pixel as active.
    //First Get the deadPixel List
    G4bool thisPixelisAlive=true;
    std::vector<G4int>DeadPxL = getDeadPixelList(IHpdNum,Isector);
    for(size_t ideadP=0; ideadP<DeadPxL.size(); ideadP++){
      if(ipixel == DeadPxL[ideadP] )thisPixelisAlive=false;
    }
    PixelIsAlive= thisPixelisAlive;
    if(thisPixelisAlive){
      RichTbRODummySD * DummySensi = new RichTbRODummySD;
    
      HpdSiPxLog->SetSensitiveDetector(DummySensi);
    }
    RichTbHpdSiPixelLVol=HpdSiPxLog;
    RichTbHpdSiPixelPVol=HpdSiPxPhys;


}






