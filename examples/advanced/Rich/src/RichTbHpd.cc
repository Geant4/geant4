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
// RichTbHpd.cc for Rich of LHCb
// History:
// Created: Sajan Easo (Sajan.Easo@cern.ch)
// Revision and changes: Patricia Mendez (Patricia.Mendez@cern.ch)
/////////////////////////////////////////////////////////////////////////////
#include <iostream>
#include "globals.hh"
#include "G4Box.hh"
#include "G4Trd.hh"
#include "G4Cons.hh"
#include "G4Tubs.hh"
#include "G4Sphere.hh"
#include "G4RotationMatrix.hh"
#include "G4ThreeVector.hh"
#include "G4Transform3D.hh"
#include "G4LogicalVolume.hh"
#include "G4VPhysicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4UnionSolid.hh"
#include "RichTbGeometryParameters.hh"
#include "G4LogicalBorderSurface.hh"
#include "G4LogicalSkinSurface.hh"
#include "G4OpticalSurface.hh"
#include "G4SDManager.hh"
#include "RichTbSD.hh"
#include "RichTbHpd.hh"

RichTbHpd::RichTbHpd() {; }
RichTbHpd::~RichTbHpd() { ; }

RichTbHpd::RichTbHpd(RichTbMaterial* RMaterial, 
                      G4VPhysicalVolume* MotherOfHpd, G4int ihpd,
		     G4int NumberofSectorsInHpd,
                      G4bool ConstructTrackingGeomSwitch){

  HpdNumber=ihpd;
  NumSectInHpd= NumberofSectorsInHpd;
  ConstructTrackGeomSwitch=ConstructTrackingGeomSwitch;

   G4Tubs* HpdMaster=
         new G4Tubs("HPDMaster",HpdMasterInnerRad,
                    HpdMasterRad,HpdMasterHalfZ,
                     HpdMasterStartPhi,HpdMasterEndPhi);

  //Overall rotations of the HPD 
    G4RotationMatrix HpdMasterRotationX,HpdMasterRotationY;
    G4RotationMatrix HpdMasterRotationZ;
    HpdMasterRotationZ.rotateZ(HpdMasterRotZ[ihpd]);
    G4ThreeVector HpdMasterTrsl(HpdMasterPosX[ihpd],HpdMasterPosY[ihpd],
                                HpdMasterPosZ[ihpd]);
    G4Transform3D HpdMasterTransform(HpdMasterRotationZ,
                                     HpdMasterTrsl);

    //Now for the placements of the HpdMaster

    // use vacuum . the photoelectrons travel in vacuum. hence
    // just for simplicity the the hpdmaster is set as vacuum.
     G4LogicalVolume* HpdMasterLog=
      new G4LogicalVolume(HpdMaster,RMaterial->getVacuum(),
                        "HpdMaster",0,0,0);
    G4VPhysicalVolume* HpdMasterPhys=
           new G4PVPlacement(HpdMasterTransform,"HpdMaster", 
                HpdMasterLog, MotherOfHpd ,false,HpdNumber);

    //Now create the Silicon Detector sectors inside the HPD.

     for (G4int isec=0; isec<NumberOfSiDetSectors; isec++){   
   
    richSiliconDetectorSect[isec] = new RichTbSiDet(RMaterial,HpdMasterPhys,
                             ConstructTrackGeomSwitch,ihpd,isec);

     }
    //Now for the rest of the components of the HPD

   if(ConstructTrackingGeomSwitch) {

    G4Tubs* HpdEnvelopeLargeTube=
      new G4Tubs("HpdEnvelopeLargeTube",HpdEnvelopeLargeTubeInnerRad,
                 HpdEnvelopeLargeTubeOuterRad,HpdEnvelopeLargeTubeHalfZ,
                 HpdEnvelopeLargeTubeStartPhi, HpdEnvelopeLargeTubeEndPhi);

    G4Cons* HpdEnvelopeCone =
      new G4Cons("HpdEnvelopeCone",HpdEnvelopeConeInnerR1,
		 HpdEnvelopeConeOuterR1,HpdEnvelopeConeInnerR2,
                 HpdEnvelopeConeOuterR2,HpdEnvelopeConeHalfZ,
                 HpdEnvelopeConeStartPhi,HpdEnvelopeConeEndPhi);

    G4Tubs* HpdEnvelopeSmallTube=
      new G4Tubs("HpdEnvelopeSmallTube",HpdEnvelopeSmallTubeInnerRad,
                 HpdEnvelopeSmallTubeOuterRad,HpdEnvelopeSmallTubeHalfZ,
                 HpdEnvelopeSmallTubeStartPhi, HpdEnvelopeSmallTubeEndPhi);

    G4Tubs* HpdEnvelopeEndCap=
      new G4Tubs("HpdEnvelopeEndcap",HpdEnvelopeEndCapInnerRad,
                 HpdEnvelopeEndCapRad,
		 HpdEnvelopeEndCapHalfZ,HpdEnvelopeEndCapStartPhi,
                 HpdEnvelopeEndCapEndPhi);

    G4Sphere* HpdQuartzWindowSphereseg =
      new G4Sphere("HpdQuartzWindow", HpdQuarzWindowRInner,
	    HpdQuarzWindowROuter,HpdQuartzWStartPhi,
            HpdQuartzWPhiSize,HpdQuartzWStartTheta,
            HpdQuartzWThetaSize);

    G4Sphere* HpdPhCathodeSphereseg =
      new G4Sphere("HpdPhCathodeWindow", HpdPhCathodeRInner,
	    HpdPhCathodeROuter,HpdPhCathodeStartPhi,
            HpdPhCathodePhiSize,HpdPhCathodeStartTheta,
            HpdPhCathodeThetaSize);

    //rotations and translations
   //Define Rotations for the HPD components 
  //                 even if they have no real rotations.

    G4RotationMatrix HpdEnvelopeLargeTubeRotation,
                     HpdEnvelopeConeRotation,
                     HpdEnvelopeSmallTubeRotation;
    G4RotationMatrix HpdEnvelopeEndCapRotation,HpdQuartzWRotation;
    G4RotationMatrix HpdPhCathodeRotation;


    G4ThreeVector HpdEnvelopeLargeTubeTrsl(HpdEnvelopeLargeTubePosX,
					HpdEnvelopeLargeTubePosY,
					HpdEnvelopeLargeTubePosZ);
    G4Transform3D HpdEnvelopeLargeTubeTransform(HpdEnvelopeLargeTubeRotation,
				               HpdEnvelopeLargeTubeTrsl);

    G4ThreeVector HpdEnvelopeConeTrsl(HpdEnvelopeConeShiftX,
				      HpdEnvelopeConeShiftY,
				      HpdEnvelopeConeShiftZ);
    G4Transform3D HpdEnvelopeConeTransform(HpdEnvelopeConeRotation,
				           HpdEnvelopeConeTrsl);


    G4ThreeVector HpdEnvelopeSmallTubeTrsl(HpdEnvelopeSmallTubeShiftX,
					HpdEnvelopeSmallTubeShiftY,
					HpdEnvelopeSmallTubeShiftZ);
    G4Transform3D HpdEnvelopeSmallTubeTransform(HpdEnvelopeSmallTubeRotation,
				               HpdEnvelopeSmallTubeTrsl);

    G4ThreeVector HpdEnvelopeEndCapTrsl(HpdEnvelopeEndCapShiftX,
				     HpdEnvelopeEndCapShiftY,
				     HpdEnvelopeEndCapShiftZ);
    G4Transform3D HpdEnvelopeEndCapTransform(HpdEnvelopeEndCapRotation,
				               HpdEnvelopeEndCapTrsl);

    G4ThreeVector HpdQuartzWTrsl(HpdQuartzWPosX,
				 HpdQuartzWPosY,
				 HpdQuartzWPosZ);
    G4Transform3D HpdQuartzWTransform(HpdQuartzWRotation,
				               HpdQuartzWTrsl);
    G4ThreeVector HpdPhCathodeTrsl(HpdPhCathodePosX,
				   HpdPhCathodePosY,
				   HpdPhCathodePosZ);
    G4Transform3D HpdPhCathodeTransform(HpdPhCathodeRotation,
				               HpdPhCathodeTrsl);

    //Now for the placements
    G4UnionSolid* HpdEnvelopeTubeLC =
      new G4UnionSolid("HpdEnvelopeTubeLC", HpdEnvelopeLargeTube, 
                       HpdEnvelopeCone,
		       HpdEnvelopeConeTransform);
     
//    G4LogicalVolume*  HpdEnvelopeTubeLCLog=
      new G4LogicalVolume(HpdEnvelopeTubeLC,RMaterial->getHpdTubeMaterial(),
			  "HpdEnvelopeTubeLC",0,0,0);

    G4UnionSolid* HpdEnvelopeTubeLCS =
      new G4UnionSolid("HpdEnvelopeTubeLCS", HpdEnvelopeTubeLC, 
                       HpdEnvelopeSmallTube,
		       HpdEnvelopeSmallTubeTransform);
     
//    G4LogicalVolume*  HpdEnvelopeTubeLCSLog=
      new G4LogicalVolume(HpdEnvelopeTubeLCS,RMaterial->getHpdTubeMaterial(),
			  "HpdEnvelopeTubeLCS",0,0,0);

    G4UnionSolid* HpdEnvelopeTubeGrand =
      new G4UnionSolid("HpdEnvelopeTubeGrand", HpdEnvelopeTubeLCS, 
                      HpdEnvelopeEndCap,HpdEnvelopeEndCapTransform);


    G4LogicalVolume*  HpdEnvelopeTubeGrandLog=
      new G4LogicalVolume(HpdEnvelopeTubeGrand,RMaterial->getHpdTubeMaterial(),
			  "HpdEnvelopeTubeGrand",0,0,0);

    //
    //
    G4LogicalVolume* HpdQuartzWindowLog=
      new G4LogicalVolume( HpdQuartzWindowSphereseg,
			   RMaterial->getPadHpdQuartzWMaterial(),
                           "HpdQuartzWindow",0,0,0);

   G4LogicalVolume* HpdPhCathodeLog=
      new G4LogicalVolume( HpdPhCathodeSphereseg,
			   RMaterial->getPadHpdPhCathodeMaterial(),
                           "HpdPhCathode",0,0,0);
   //
   //
    G4VPhysicalVolume* HpdEnvelopeTubeGrandPhys=
           new G4PVPlacement(HpdEnvelopeLargeTubeTransform,
                "HpdEnvelopeTubeGrand", 
                HpdEnvelopeTubeGrandLog,HpdMasterPhys,false,HpdNumber);

    G4VPhysicalVolume* HpdQuartzWindowPhys=
           new G4PVPlacement(HpdQuartzWTransform,"HpdQuartzWindow", 
                HpdQuartzWindowLog,HpdMasterPhys,false,HpdNumber);

    G4VPhysicalVolume* HpdPhCathodePhys=
           new G4PVPlacement(HpdPhCathodeTransform,"HpdPhCathode", 
                HpdPhCathodeLog,HpdMasterPhys,false,HpdNumber);


//   G4LogicalBorderSurface* HpdQuartzWSurface =
     new G4LogicalBorderSurface("HpdQuartzWSurface",HpdMasterPhys,
                                 HpdQuartzWindowPhys,
			  RMaterial->getOpticalHpdTQuartzWSurface());

//   G4LogicalBorderSurface* HpdQPhCathodeSurface =
     new G4LogicalBorderSurface("HpdQPhCathodeSurface",HpdQuartzWindowPhys,
                                 HpdPhCathodePhys,
		 RMaterial->getOpticalHpdQuartzWPhCathodeSurface());

//   G4LogicalBorderSurface* HpdPhCathodeQSurface =
     new G4LogicalBorderSurface("HpdQPhCathodeSurface", HpdPhCathodePhys,
                                HpdQuartzWindowPhys,
		 RMaterial->getOpticalHpdQuartzWPhCathodeSurface());

   //change made on 26-9-01 SE to make photocathode transparent.
   //    G4LogicalSkinSurface* HpdPhCathodeSkinSurface =
   //    new G4LogicalSkinSurface("HpdPhCathodeSkinSurface",
   //			       HpdPhCathodeLog,
   //		 RMaterial->getOpticalPhCathodeSkinSurface());

//    G4LogicalBorderSurface* HpdPhCathodeMasterSurface =
      new G4LogicalBorderSurface("HpdPhCathodeMasterSurface",
				 HpdPhCathodePhys,HpdMasterPhys,   
		 RMaterial->getOpticalPhCathodeBorderSurface());

//    G4LogicalBorderSurface* HpdMasterPhCathodeSurface =
      new G4LogicalBorderSurface("HpdMasterPhCathodeSurface",
			HpdMasterPhys,HpdPhCathodePhys,   
		 RMaterial->getOpticalPhCathodeBorderSurface());

//   G4LogicalBorderSurface* HpdEnvelopeSurface =
     new G4LogicalBorderSurface("HpdEnvelopeSurface",HpdMasterPhys,
                          HpdEnvelopeTubeGrandPhys,
				RMaterial->getOpticalHpdMetalSurface());
    RichTbHpdEnvelopeTubeGrandLVol=HpdEnvelopeTubeGrandLog;
    RichTbHpdEnvelopeTubeGrandPVol=HpdEnvelopeTubeGrandPhys;
    RichTbHpdQuartzWLVol=HpdQuartzWindowLog;
    RichTbHpdQuartzWPVol=HpdQuartzWindowPhys;
    RichTbPhCathodeLVol=HpdPhCathodeLog;
    RichTbPhCathodePVol=HpdPhCathodePhys;

   }
    RichTbHpdLVol=HpdMasterLog;
    RichTbHpdPVol=HpdMasterPhys;
       

}







