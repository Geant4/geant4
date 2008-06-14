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
// RichTbGeometryParameters.hh for Rich of LHCb
// History:
// Created: Sajan Easo (Sajan.Easo@cern.ch)
// Revision and changes: Patricia Mendez (Patricia.Mendez@cern.ch)
/////////////////////////////////////////////////////////////////////////////
#ifndef RichTbGeometryParameters_h
#define RichTbGeometryParameters_h 1

#include "globals.hh"
#include "AerogelTypeSpec.hh"
#include <cmath>
extern void InitializeRichTbGeometry();
extern G4double GetCurAerogelLength(G4int);
//
//
static const G4double sqroot3=std::pow(3.0,0.5);

//Size of the LHCb Rich Testbeam Hall. 
static const G4double ExpHallHalfX=4000.0*mm;
static const G4double ExpHallHalfY=4000.0*mm;
static const G4double ExpHallHalfZ=8000.0*mm;
// The Hall is kept at the orgin of the coord system.
// The coord system has +z along the beam direction and +y
// going upwards.
//
//Now for the Vessel
static const G4double VesselInnerRad=0.0*mm;
static const G4double VesselOuterRad=275.0*mm;
static const G4double VesselHalfZ=540.0*mm;
static const G4double VesselStartPhi=0.0*rad;
static const G4double VesselDelPhi=twopi*rad;
static const G4double VesselPosX=0.0*mm;
static const G4double VesselPosY=0.0*mm;
static const G4double VesselPosZ=VesselHalfZ;

//Now for the box containing the  aerogel sample.
static const G4double RadFrameHalfX=55.0*mm;
static const G4double RadFrameHalfY=55.0*mm;
static const G4double RadFrameHalfZ=55.0*mm;
static const G4double RadFramePosX=0.0*mm;
static const G4double RadFramePosY=0.0*mm;
//The following are the Z positons of start of the radiator frame
// and aerogel tiles.
static const G4double RadFrameGenStartZ=185.0*mm;
static const G4double AgelTileGenStartZ=190.0*mm;
static const G4double RadFramePosZ=
            RadFrameGenStartZ+RadFrameHalfZ-VesselHalfZ;
// Now for the aerogel blocks. They are made to start at 190mm from the
// upstream edge of the box.
// The following are dimensioned to MaxNumberOfAerogelTiles which
// from the AerogelTypeSpec.hh file
//
//Now for the window of the radframe at its upstream and downstream
// ends.
//
static const G4double RadHoldUpHalfX=50.0*mm;
static const G4double RadHoldUpHalfY=50.0*mm;
static const G4double RadHoldUpHalfZ=1.0*mm;
static const G4double RadHoldUpPosX=0.0*mm;
static const G4double RadHoldUpPosY=0.0*mm;
static const G4double RadHoldUpPosZ=-RadFrameHalfZ+ RadHoldUpHalfZ;
static const G4double RadWinUpOuterRad=40.0*mm;
static const G4double RadWinUpInnerRad=0.0*mm;
static const G4double RadWinUpHalfZ= RadHoldUpHalfZ+5.0*mm;
static const G4double RadWinUpStartPhi=0.0*rad;
static const G4double RadWinUpDelPhi=twopi*rad;
static const G4double RadWinUpShiftX=0.0*mm;
static const G4double RadWinUpShiftY=0.0*mm;
static const G4double RadWinUpShiftZ=0.0*mm;
//
static const G4double RadHoldDnHalfX=50.0*mm;
static const G4double RadHoldDnHalfY=50.0*mm;
static const G4double RadHoldDnHalfZ=1.0*mm;
static const G4double RadHoldDnPosX=0.0*mm;
static const G4double RadHoldDnPosY=0.0*mm;
static const G4double RadHoldDnShiftZ=15.0*mm;
static const G4double RadHoldDnPosZ=RadFrameHalfZ- 
                      RadHoldDnShiftZ- RadHoldDnHalfZ;
static const G4double RadWinDnOuterRad=40.0*mm;
static const G4double RadWinDnInnerRad=0.0*mm;
static const G4double RadWinDnHalfZ= RadHoldDnHalfZ+5.0*mm;
static const G4double RadWinDnStartPhi=0.0*rad;
static const G4double RadWinDnDelPhi=twopi*rad;
static const G4double RadWinDnShiftX=0.0*mm;
static const G4double RadWinDnShiftY=0.0*mm;
static const G4double RadWinDnShiftZ=0.0*mm;
// the following three variables are for each of the aerogel types. 
// In the G4example only 1 type is simply repeated 5 times.
static const G4double AgelHalfX[]={35.0*mm,35.0*mm,35.0*mm,35.0*mm,35.0*mm};
static const G4double AgelHalfY[]={40.0*mm,40.0*mm,40.0*mm,40.0*mm,40.0*mm};
static const G4double AgelHalfZ[]={20.0*mm,20.0*mm,20.0*mm,20.0*mm,20.0*mm};
// The following 2 variables are for each tile.
// for now no shifts forseen in the XY direction.
static const G4double AgelPosX[]={0.0*mm,0.0*mm,0.0*mm,0.0*mm,0.0*mm};
static const G4double AgelPosY[]={0.0*mm,0.0*mm,0.0*mm,0.0*mm,0.0*mm};
static const G4double AgelStartZ= 
                AgelTileGenStartZ-RadFrameGenStartZ-RadFrameHalfZ;
static const G4double TotalAgelThickness=80.0*mm;
static const G4double AgelTileGapZ=1.0*mm;
static const G4double AgelEndZ= AgelStartZ+TotalAgelThickness+AgelTileGapZ;


//Now for the wraps above and below the aerogel tiles.
static const G4double AgelWrapTopHalfX[]=
               {35.0*mm,35.0*mm,35.0*mm,35.0*mm,35.0*mm};
static const G4double AgelWrapTopHalfY[]=
               {1.0*mm,1.0*mm,1.0*mm,1.0*mm,1.0*mm};
static const G4double AgelWrapTopHalfZ[]=
               {20.0*mm,20.0*mm,20.0*mm,20.0*mm,20.0*mm};

static const G4double AgelWrapTopPosX[]=
               {0.0*mm,0.0*mm,0.0*mm,0.0*mm,0.0*mm};

static const G4double AgelWrapBotHalfX[]=
               {35.0*mm,35.0*mm,35.0*mm,35.0*mm,35.0*mm};
static const G4double AgelWrapBotHalfY[]=
               {1.0*mm,1.0*mm,1.0*mm,1.0*mm,1.0*mm};
static const G4double AgelWrapBotHalfZ[]=
               {20.0*mm,20.0*mm,20.0*mm,20.0*mm,20.0*mm};

static const G4double AgelWrapBotPosX[]=
               {0.0*mm,0.0*mm,0.0*mm,0.0*mm,0.0*mm};
// in the G4Example only 1 type of filter is used.
static const G4double FilterHalfX=53.0*mm;
static const G4double FilterHalfY=53.0*mm;
static const G4double GlassD263HalfZ=0.15*mm;
static const G4double FilterHalfZArray[]={GlassD263HalfZ, GlassD263HalfZ,
       GlassD263HalfZ, GlassD263HalfZ,  GlassD263HalfZ,   GlassD263HalfZ };
// have a nominal value for the filter thickness.
static const G4double FilterHalfZNominal= GlassD263HalfZ;
static const G4double FilterPosX=0.0*mm;
static const G4double FilterPosY=0.0*mm;
//gap between aerogel and Filter in Z.
static const G4double FilterAgelGapZ=2.5*mm;
//nominal value for the Filter position
static const G4double FilterPosZNominal=
                 AgelEndZ+FilterAgelGapZ+ FilterHalfZNominal; 
//Now for the mirror
//
static const G4double MirrorRInner=1185.0*mm;
static const G4double MirrorROuter=1191.0*mm;
static const G4double MirrorHorizontalChord=350.0*mm;
static const G4double MirrorVerticalChord=290.0*mm;
// the following is 600+117+6 mm.
// here 121 is the hpd Q window outer Z and 6 is the
// thickness of the mirror 
// static const G4double MirrorShiftFromEnd=707.0*mm;
static const G4double MirrorShiftFromEnd=723.0*mm;
static const G4double MirrorNominalPosZ=MirrorShiftFromEnd-VesselHalfZ
                                 -MirrorRInner;
static const G4double MirrorNominalRotX=0.0*rad;
static const G4double MirrorNominalRotY=0.0*rad;

//Now for each of the HPDs.
// First the size of each part.
static const   G4double HpdMasterRad=64.00*mm;
//static const   G4double HpdMasterHalfZ=55.0*mm;
static const   G4double HpdMasterHalfZ=60.0*mm;
static const   G4double HpdEnvelopeLargeTubeInnerRad=58.5*mm;
static const   G4double HpdEnvelopeLargeTubeOuterRad=63.5*mm;
//static const   G4double HpdEnvelopeLargeTubeHalfZ=8.73*mm;
static const   G4double HpdEnvelopeLargeTubeHalfZ=11.0*mm;
static const   G4double HpdEnvelopeConeHalfZ=20.0*mm;
static const   G4double HpdEnvelopeSmallTubeHalfZ=10.0*mm;
static const   G4double HpdEnvelopeSmallTubeInnerRad=35.0*mm;
static const   G4double HpdEnvelopeSmallTubeOuterRad=40.0*mm;
static const   G4double HpdEnvelopeConeOuterR2=HpdEnvelopeLargeTubeOuterRad;
static const   G4double HpdEnvelopeConeInnerR2=HpdEnvelopeLargeTubeInnerRad;
static const   G4double HpdEnvelopeConeOuterR1=HpdEnvelopeSmallTubeOuterRad;
static const   G4double HpdEnvelopeConeInnerR1=HpdEnvelopeSmallTubeInnerRad;
//
static const   G4double HpdEnvelopeEndCapRad=HpdEnvelopeSmallTubeOuterRad;
static const   G4double HpdEnvelopeEndCapHalfZ=2.5*mm;
static const   G4double HpdQuartzWindowThickness=4.0*mm;
static const   G4double HpdQuarzWindowROuter=100.0*mm;
static const   G4double HpdQuarzWindowRInner= HpdQuarzWindowROuter
                         -HpdQuartzWindowThickness;
static const G4double PhotoCathodeThickness=0.00004*mm;
static const   G4double HpdPhCathodeROuter= HpdQuarzWindowRInner;
static const   G4double HpdPhCathodeRInner= 
             HpdQuarzWindowRInner - PhotoCathodeThickness;
// The following obtained by arcsin(63.5/100)
static const   G4double HpdQuartzWThetaSize=0.6880*rad;
// The following obtained by arcsin(57/(100-4))
static const   G4double HpdPhCathodeThetaSize=0.6357*rad;

// For defining cylinders and spherical segment..
static const   G4double HpdMasterInnerRad=0.0*mm;
static const   G4double HpdMasterStartPhi=0.0*rad;
static const   G4double HpdMasterEndPhi=twopi*rad;

static const   G4double HpdEnvelopeLargeTubeStartPhi=0.0*rad;
static const   G4double HpdEnvelopeLargeTubeEndPhi=twopi*rad;

static const   G4double HpdEnvelopeConeStartPhi=0.0*rad;
static const   G4double HpdEnvelopeConeEndPhi=twopi*rad;


static const   G4double HpdEnvelopeSmallTubeStartPhi=0.0*rad;
static const   G4double HpdEnvelopeSmallTubeEndPhi=twopi*rad;

static const   G4double HpdEnvelopeEndCapInnerRad=0.0*mm;
static const   G4double HpdEnvelopeEndCapStartPhi=0.0*rad;
static const   G4double HpdEnvelopeEndCapEndPhi=twopi*rad;

static const   G4double HpdQuartzWPhiSize=twopi*rad;
static const   G4double HpdQuartzWStartTheta=0.0*rad;
static const   G4double HpdQuartzWStartPhi=0.0*rad;
static const   G4double HpdPhCathodePhiSize=twopi*rad;
static const   G4double HpdPhCathodeStartTheta=0.0*rad;
static const   G4double HpdPhCathodeStartPhi=0.0*rad;
// Z locations of the various parts.
static const   G4double HpdQuartzPartFromEndZ=4.0*mm;
//static const   G4double HpdEnvelopePartFromEndZ=18.0*mm;
// The following obtained by requiring 100-std::sqrt(100*100-63.5*63.5)=22.75 
// for the quartz region in Z. The total is 22.75+4+4=30.75.  
static const   G4double HpdEnvelopePartFromEndZ=30.75*mm;
static const G4double HpdPhotoCathodeSiZdist=100.0*mm;


//Silicon detector inside the HPD
static const G4int NumberOfSiDetSectors=16;
// Now for each one of the sector.
// In the following the 0.0001 is to avoid graphics from crashing
// for zero length of a triangle made from a trapozoid.

static const G4double SiSectAngSize=(22.5*pi/180)*rad;
static const G4double SiSectAngHalfSize=SiSectAngSize/2.0;
static const G4double SiSectHeight=25.0*mm;
static const G4double SiSectSide= SiSectHeight/std::cos(SiSectAngHalfSize);
static const G4double SiSectHalfMoonGap=0.0*mm;
static const G4double SiSectTrapHalfY1=0.15*mm;  // halfthickness of Si.
static const G4double SiSectTrapHalfY2= SiSectTrapHalfY1;
static const G4double SiSectTrapHalfX1=0.0001*mm;
static const G4double SiSectTrapHalfX2 = SiSectHeight*std::tan(SiSectAngHalfSize);
static const G4double SiSectTrapHalfZ  = SiSectHeight/2.0 ;
static const G4double SiSectRotX=(90*pi/180)*rad;
static const G4double SiSectAngStart=(pi/2.0)*rad-SiSectAngSize/2.0;
// SiSectPosX and SiSectPosY are calculated in the cc file.
static const G4double SiSectPosZ=
         HpdMasterHalfZ-HpdQuartzPartFromEndZ- 
         HpdQuartzWindowThickness-HpdPhotoCathodeSiZdist;
//
//
//Now for the coating on the Silicon surface.
static const G4double SiSectCoatingAngSize=SiSectAngSize;
static const G4double SiSectCoatingAngHalfSize=SiSectAngSize/2.0;
static const G4double SiSectCoatingSide=SiSectSide;
static const G4double SiSectCoatingHalfMoonGap=SiSectHalfMoonGap;
static const G4double SiSectCoatingTrapHalfY1=0.05*mm; // SiCoating HalfThick.
static const G4double SiSectCoatingTrapHalfY2= SiSectCoatingTrapHalfY1;
static const G4double SiSectCoatingTrapHalfX1= SiSectTrapHalfX1;
static const G4double SiSectCoatingTrapHalfX2= SiSectTrapHalfX2;
static const G4double SiSectCoatingTrapHalfZ= SiSectTrapHalfZ;
static const G4double SiSectCoatingRotX=SiSectRotX;
static const G4double  SiSectCoatingAngStart=SiSectAngStart;
// SiSectCoatingPosX and SiSectCoatingPosY are calculated in the cc file.
static const G4double SiSectCoatingPosZ =
      SiSectPosZ + SiSectTrapHalfY1+SiSectCoatingTrapHalfY1;
//
//
static const G4double XsizePix=1.0*mm;
static const G4double YsizePix=1.0*mm;
static const G4double YsizeBigPix=2.0*mm;
// the number of the big pixel in the pixel map on the web=59, hence
// in the c++ array 58 since the array start from 0.
static const G4int BigPixelNum=58;
//Silicon pixels inside the Silicon detector.
static const G4double SiPixelHalfY= SiSectTrapHalfY1;
static const G4double SiPixelHalfX= XsizePix/2.0;
static const G4double SiPixelHalfZ= YsizePix/2.0;
// now for the big pixel at the centre of the hpd.
static const G4double SiBigPixelHalfZ=YsizeBigPix/2.0;
//Now for their relative locations
static const G4double HpdEnvelopeLargeTubePosX=0.0*mm;
static const G4double HpdEnvelopeLargeTubePosY=0.0*mm;
static const G4double HpdEnvelopeLargeTubePosZ=
         HpdMasterHalfZ-HpdEnvelopePartFromEndZ-HpdEnvelopeLargeTubeHalfZ;

static const G4double HpdEnvelopeConeShiftX=0.0*mm;
static const G4double HpdEnvelopeConeShiftY=0.0*mm;
static const G4double HpdEnvelopeConeShiftZ=
    -(HpdEnvelopeLargeTubeHalfZ+HpdEnvelopeConeHalfZ);

static const G4double HpdEnvelopeSmallTubeShiftX=0.0*mm;
static const G4double HpdEnvelopeSmallTubeShiftY=0.0*mm;
static const G4double HpdEnvelopeSmallTubeShiftZ=
            -(HpdEnvelopeLargeTubeHalfZ+2*HpdEnvelopeConeHalfZ+
              HpdEnvelopeSmallTubeHalfZ);
static const G4double HpdEnvelopeEndCapShiftX=0.0*mm;
static const G4double HpdEnvelopeEndCapShiftY=0.0*mm;
static const G4double HpdEnvelopeEndCapShiftZ=
-(HpdEnvelopeLargeTubeHalfZ+2*HpdEnvelopeConeHalfZ+
  2*HpdEnvelopeSmallTubeHalfZ+HpdEnvelopeEndCapHalfZ);

static const   G4double HpdQuartzWPosX=0.0*mm;
static const   G4double HpdQuartzWPosY=0.0*mm;
static const   G4double HpdQuartzWPosZ= 
       HpdMasterHalfZ- HpdQuartzPartFromEndZ-HpdQuarzWindowROuter;

static const   G4double HpdPhCathodePosX=0.0*mm;
static const   G4double HpdPhCathodePosY=0.0*mm;
static const   G4double HpdPhCathodePosZ= HpdQuartzWPosZ; 


//Placement of the Hpd Si pixels in the Hpd Si Sector.
static const G4int NumberOfPadHpdSiPixels=128;
static const G4int MaxNumberOfPixRow=24;
static const G4int MaxNumberOfPixCol=10;
extern G4double PixRowNumSect[NumberOfPadHpdSiPixels];
extern G4double PixColNumSect[NumberOfPadHpdSiPixels];
extern G4bool PixelAtSectEdge[NumberOfPadHpdSiPixels];
// the following is for the rows after the central big pixel.
// the central pixel is at row 0. 
// the others start from row number start at 1.
// the following is 1 mm in the current setup.
static const G4double RowInitPointBigPixel= SiSectTrapHalfZ*2.0
                               -(MaxNumberOfPixRow-1) * YsizePix
                               -  YsizeBigPix/2.0;
// the following is 1.5 mm in the current setup.
static const G4double RowInitPoint= SiSectTrapHalfZ*2.0
                               -(MaxNumberOfPixRow-1) * YsizePix
                               - YsizePix/2.0;
//The follwing two variables are defined in the cc file for this
//include file as part of the class declared below. Since the
// Si det plane is in the XZ plane, a swap of Y->Z and Z->Y is done
// while positioning the Si Pixel.
static const G4double HpdSiPixPosZ=0.0*mm;
#ifndef RichTbPadHpdSiPixPos_h
#define RichTbPadHpdSiPixPos_h 1
class RichTbPadHpdSiPixPos{
public:
  RichTbPadHpdSiPixPos(G4int);
  virtual ~RichTbPadHpdSiPixPos();
  G4double getPadHpdSiPixPosX() {return PadHpdSiPixPosX;}
  G4double getPadHpdSiPixPosY() {return PadHpdSiPixPosY;}
  G4double getCurrentPixelnum(){return icurpixel;}  
private:
G4double PadHpdSiPixPosX;
G4double PadHpdSiPixPosY;
G4int icurpixel;
};
#endif 

//Placement of all the HPDs.

static const G4int NumberOfHpds=4;
static const G4double HpdPosRad=146.5*mm;
static const G4double HpdZfromEnd=HpdMasterHalfZ-VesselPosZ+1.0*mm;

static const  G4double HpdMasterPosX[NumberOfHpds]=
    {HpdPosRad,0.0*mm,-HpdPosRad,0.0*mm};
static const  G4double HpdMasterPosY[NumberOfHpds]=
     {0.0*mm,-HpdPosRad,0.0*mm,HpdPosRad};
static const  G4double HpdMasterPosZ[NumberOfHpds]=
     {HpdZfromEnd,HpdZfromEnd,HpdZfromEnd,HpdZfromEnd};
//rot1 version of hpd rotations.
static const G4double HpdMasterRotZ[NumberOfHpds]=
 {0.0*rad,(69.5*pi/180.0)*rad,(22.5*pi/180.0)*rad,(7.5*pi/180.0)*rad};
//
//
#endif 







