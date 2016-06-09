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
// RichTbGraphics.cc for Rich of LHCb
// History:
// Created: Sajan Easo (Sajan.Easo@cern.ch)
// Revision and changes: Patricia Mendez (Patricia.Mendez@cern.ch)
/////////////////////////////////////////////////////////////////////////////
#include <iostream>
#include "RichTbDetectorConstruction.hh"
#include "RichTbHall.hh"
#include "RichTbComponent.hh"
#include "RichTbGeometryParameters.hh"
#include "RichTbPhotoDetector.hh"
#include "RichTbHpd.hh"
#include "RichTbSiDet.hh"
#include "RichTbGraphics.hh"
#include "G4LogicalVolume.hh"
#include "G4VisAttributes.hh"
#include "G4Colour.hh"
#include "RichTbRunConfig.hh"
RichTbGraphics::RichTbGraphics(){ ; }
RichTbGraphics::RichTbGraphics(RichTbHall* RTbHall,
    RichTbComponent* RTbComponent,
    RichTbPhotoDetector* RTbPhotoDetector, 
    RichTbROGeometry* RTbROGeom,
    RichTbRunConfig* RConfig){

  //For Following variables 0 means make the volume invisible;
  //                        1 means make it visible as a solid.
  //                        2 means make it visible as a wireframe.
  G4int RichTbHall_vis= RConfig->getRichTbHall_visib(); 
  G4int RichTbEnclosure_vis = RConfig->getRichTbEnclosure_visib() ;
  G4int RichTbRadFrame_vis = RConfig->getRichTbRadFrame_visib() ;
  G4int RichTbRadUpW_vis = RConfig->getRichTbRadUpW_visib() ;
  G4int RichTbRadDnW_vis = RConfig->getRichTbRadDnW_visib() ;
  G4int RichTbAerogel_vis = RConfig->getRichTbAerogel_visib() ;
  G4int RichTbAerogelWrap_vis = RConfig->getRichTbAerogelWrap_visib() ;
  G4int RichTbAerogelWrapTop_vis = RichTbAerogelWrap_vis;
  G4int RichTbAerogelWrapBot_vis =  RichTbAerogelWrap_vis;
  G4int RichTbFilter_vis = RConfig->getRichTbFilter_visib() ;
  G4int RichTbMirror_vis =   RConfig->getRichTbMirror_visib();
  G4int RichTbHpdMaster_vis=  RConfig->getRichTbHpdMaster_visib() ;
  G4int RichTbHpdEnvelopeTube_vis=  RConfig->getRichTbHpdEnvelopeTube_visib();
  G4int RichTbHpdQuartzW_vis= RConfig->getRichTbHpdQuartzW_visib() ;
  G4int RichTbHpdPhCathode_vis= RConfig->getRichTbHpdPhCathode_visib() ;
  G4int RichTbHpdSiDet_vis= RConfig->getRichTbHpdSiDet_visib() ;
  G4int RichTbHpdSectCoat_vis =RConfig->getRichTbHpdSectCoat_visib() ;
  G4int RichTbHpdSiPx_vis= RConfig->getRichTbHpdSiPx_visib() ;
  
  //
  //Now for the RichTbHall

  G4LogicalVolume* RichTbHall_LV = RTbHall->getRichTbHallLogicalVolume();
  if(RichTbHall_vis == 0 ) {
     RichTbHall_LV->SetVisAttributes(G4VisAttributes::Invisible);
  }else {
     G4VisAttributes* RichTbHall_logVisAtt 
     = new G4VisAttributes(G4Colour(0.,1.,1.));
     if(RichTbHall_vis == 2 ){
        RichTbHall_logVisAtt->SetForceWireframe(true);
     }
     RichTbHall_logVisAtt->SetVisibility(true);
     RichTbHall_LV->SetVisAttributes(RichTbHall_logVisAtt);
  //
  }
  //
  //
  //Now for the RichTbEnclosure

  G4LogicalVolume* RichTbEnclosure_LV = 
           RTbComponent->getEnclosureLogicalVolume();

  if(RichTbEnclosure_vis == 0 ){
     RichTbEnclosure_LV->SetVisAttributes (G4VisAttributes::Invisible);
  }else {

     G4VisAttributes* RichTbEnclosure_logVisAtt
      = new G4VisAttributes(G4Colour(1.,1.,0.));
     if( RichTbEnclosure_vis == 2) {
        RichTbEnclosure_logVisAtt->SetForceWireframe(true);
     }else {
        RichTbEnclosure_logVisAtt->SetForceSolid(true);
     }
     RichTbEnclosure_logVisAtt->SetVisibility(true);
     RichTbEnclosure_LV->
               SetVisAttributes(RichTbEnclosure_logVisAtt);
  } 
  //
  //Now for the Radiator frame

  G4LogicalVolume* RichTbRadFrame_LV = 
           RTbComponent->getRadFrameLogicalVolume();

  if(RichTbRadFrame_vis == 0 ){
     RichTbRadFrame_LV->SetVisAttributes (G4VisAttributes::Invisible);
  }else {

     G4VisAttributes* RichTbRadFrame_logVisAtt
      = new G4VisAttributes(G4Colour(1.,1.,0.));
     if( RichTbRadFrame_vis == 2) {
        RichTbRadFrame_logVisAtt->SetForceWireframe(true);
     }else {
        RichTbRadFrame_logVisAtt->SetForceSolid(true);
     }
     RichTbRadFrame_logVisAtt->SetVisibility(true);
     RichTbRadFrame_LV->
               SetVisAttributes(RichTbRadFrame_logVisAtt);
  } 
  //
  //Now for the upstream window of the radiator frame

  G4LogicalVolume* RichTbRadUpW_LV = 
           RTbComponent->getRadUpWLogicalVolume();

  if(RichTbRadUpW_vis == 0 ){
     RichTbRadUpW_LV->SetVisAttributes (G4VisAttributes::Invisible);
  }else {

     G4VisAttributes* RichTbRadUpW_logVisAtt
      = new G4VisAttributes(G4Colour(0.0,0.3,0.7));
     if( RichTbRadUpW_vis == 2) {
        RichTbRadUpW_logVisAtt->SetForceWireframe(true);
     }else {
        RichTbRadUpW_logVisAtt->SetForceSolid(true);
     }
     RichTbRadUpW_logVisAtt->SetVisibility(true);
     RichTbRadUpW_LV->
               SetVisAttributes(RichTbRadUpW_logVisAtt);
  } 
  //
  // Now for the downstream window of the radiator frame
  G4LogicalVolume* RichTbRadDnW_LV = 
           RTbComponent->getRadDnWLogicalVolume();

  if(RichTbRadDnW_vis == 0 ){
     RichTbRadDnW_LV->SetVisAttributes (G4VisAttributes::Invisible);
  }else {

     G4VisAttributes* RichTbRadDnW_logVisAtt
      = new G4VisAttributes(G4Colour(0.0,0.2,0.8));
     if( RichTbRadDnW_vis == 2) {
        RichTbRadDnW_logVisAtt->SetForceWireframe(true);
     }else {
        RichTbRadDnW_logVisAtt->SetForceSolid(true);
     }
     RichTbRadDnW_logVisAtt->SetVisibility(true);
     RichTbRadDnW_LV->
               SetVisAttributes(RichTbRadDnW_logVisAtt);
  } 
  //
  //
  //Now for the Aerogel

  G4int NumTiles=RConfig->GetNumberOfAerogelTiles();
  for (G4int ng=0; ng< NumTiles; ng++ ) {
  G4LogicalVolume* RichTbAerogel_LV = 
           RTbComponent->getAgelLogicalVolume(ng);

  if(RichTbAerogel_vis == 0 ){
     RichTbAerogel_LV->SetVisAttributes (G4VisAttributes::Invisible);
  }else {

    G4VisAttributes* RichTbAerogel_logVisAtt;
    if(ng == 0 || ng == 2 || ng || 4 ) {
     RichTbAerogel_logVisAtt = new G4VisAttributes(G4Colour(1.,0.2,0.8));
    }else {

    RichTbAerogel_logVisAtt = new G4VisAttributes(G4Colour(0.7,0.6,0.9));

    } 

     if( RichTbAerogel_vis == 2) {
        RichTbAerogel_logVisAtt->SetForceWireframe(true);
     }else {
        RichTbAerogel_logVisAtt->SetForceSolid(true);
     }
     RichTbAerogel_logVisAtt->SetVisibility(true);
     RichTbAerogel_LV->
               SetVisAttributes(RichTbAerogel_logVisAtt);
  }

  //Now for the top Wrap of the Aerogel tile.

  if( RTbComponent->getAgelWrapTopLogicalVolume(ng)) {
  G4LogicalVolume* RichTbAerogelWrapTop_LV = 
           RTbComponent->getAgelWrapTopLogicalVolume(ng);

  if(RichTbAerogelWrapTop_vis == 0 ){
     RichTbAerogelWrapTop_LV->SetVisAttributes (G4VisAttributes::Invisible);
  }else {

    G4VisAttributes* RichTbAerogelWrapTop_logVisAtt;
    if(ng == 0 || ng == 2 || ng || 4 ) {
     RichTbAerogelWrapTop_logVisAtt = 
              new G4VisAttributes(G4Colour(0.2,0.2,0.6));
    }else {

    RichTbAerogelWrapTop_logVisAtt = 
          new G4VisAttributes(G4Colour(0.2,0.4,0.8));

    } 

     if( RichTbAerogelWrapTop_vis == 2) {
        RichTbAerogelWrapTop_logVisAtt->SetForceWireframe(true);
     }else {
        RichTbAerogelWrapTop_logVisAtt->SetForceSolid(true);
     }
     RichTbAerogelWrapTop_logVisAtt->SetVisibility(true);
     RichTbAerogelWrapTop_LV->
               SetVisAttributes(RichTbAerogelWrapTop_logVisAtt);
  }
  }
  //Now for the Bottom Wrap of the Aerogel tile.

  if(RTbComponent->getAgelWrapBotLogicalVolume(ng)){
  G4LogicalVolume* RichTbAerogelWrapBot_LV = 
           RTbComponent->getAgelWrapBotLogicalVolume(ng);

  if(RichTbAerogelWrapBot_vis == 0 ){
     RichTbAerogelWrapBot_LV->SetVisAttributes (G4VisAttributes::Invisible);
  }else {

    G4VisAttributes* RichTbAerogelWrapBot_logVisAtt;
    if(ng == 0 || ng == 2 || ng || 4 ) {
     RichTbAerogelWrapBot_logVisAtt = 
              new G4VisAttributes(G4Colour(0.2,0.2,0.6));
    }else {

    RichTbAerogelWrapBot_logVisAtt = 
          new G4VisAttributes(G4Colour(0.2,0.4,0.8));

    } 

     if( RichTbAerogelWrapBot_vis == 2) {
        RichTbAerogelWrapBot_logVisAtt->SetForceWireframe(true);
     }else {
        RichTbAerogelWrapBot_logVisAtt->SetForceSolid(true);
     }
     RichTbAerogelWrapBot_logVisAtt->SetVisibility(true);
     RichTbAerogelWrapBot_LV->
               SetVisAttributes(RichTbAerogelWrapBot_logVisAtt);
  }

  }
 
}
  //
  //
  //Now for the Filter
  // First verify that the filter exists.
   G4int Filnum=RConfig->GetFilterTNumber() ;
   if(Filnum >=0 ) {
  
  G4LogicalVolume* RichTbFilter_LV = 
           RTbComponent->getFilterLogicalVolume();

  if(RichTbFilter_vis == 0 ){
     RichTbFilter_LV->SetVisAttributes (G4VisAttributes::Invisible);
  }else {

     G4VisAttributes* RichTbFilter_logVisAtt
      = new G4VisAttributes(G4Colour(0.3,0.9,0.4));
     if( RichTbFilter_vis == 2) {
        RichTbFilter_logVisAtt->SetForceWireframe(true);
     }else {
        RichTbFilter_logVisAtt->SetForceSolid(true);
     }
     RichTbFilter_logVisAtt->SetVisibility(true);
     RichTbFilter_LV->
               SetVisAttributes(RichTbFilter_logVisAtt);
  } 
   }
  //
  //
  //
  //Now for the Mirror
  G4LogicalVolume* RichTbMirror_LV = 
           RTbComponent->getMirrorLogicalVolume();

  if(RichTbMirror_vis == 0 ){
      RichTbMirror_LV->
           SetVisAttributes (G4VisAttributes::Invisible);
  }else {
     G4VisAttributes* RichTbMirror_logVisAtt =
      new G4VisAttributes(G4Colour(0.,0.5,0.5));
     if(RichTbMirror_vis == 2 ){
       RichTbMirror_logVisAtt->SetForceWireframe(true);
     }else{
       RichTbMirror_logVisAtt->SetForceWireframe(false);
       RichTbMirror_logVisAtt->SetForceSolid(true);
     }
     RichTbMirror_logVisAtt->SetVisibility(true);

     RichTbMirror_LV->
            SetVisAttributes(RichTbMirror_logVisAtt);
  }

  //Now for the HPDs.
  G4int numhpd= RTbPhotoDetector->getNumberOfHpds();
  for(G4int ihpd=0; ihpd<numhpd; ihpd++){
       //Now for the HPD master volume
      G4LogicalVolume* RichTbHpdMaster_LV = 
           RTbPhotoDetector->getRichHPD(ihpd)->getHpdLogicalVolume();

     if(RichTbHpdMaster_vis == 0 ){
      RichTbHpdMaster_LV->
           SetVisAttributes (G4VisAttributes::Invisible);
     }else {
 
       G4VisAttributes* RichTbHpdMaster_logVisAtt =
        new G4VisAttributes(G4Colour(0.8,0.6,0.2));
       if(RichTbHpdMaster_vis == 2 ){
         RichTbHpdMaster_logVisAtt->SetForceWireframe(true);
       }else{
         RichTbHpdMaster_logVisAtt->SetForceWireframe(false);
         RichTbHpdMaster_logVisAtt->SetForceSolid(true);
       }
       RichTbHpdMaster_logVisAtt->SetVisibility(true);

       RichTbHpdMaster_LV->
            SetVisAttributes(RichTbHpdMaster_logVisAtt);


     }
     //hpd Envelope Tube
      G4LogicalVolume* RichTbHpdEnvelopeTube_LV = 
      RTbPhotoDetector->getRichHPD(ihpd)->getHpdEnvelopeTubeGrandLogicalVolume();


     if(RichTbHpdEnvelopeTube_vis == 0 ){
      RichTbHpdEnvelopeTube_LV->
           SetVisAttributes (G4VisAttributes::Invisible);
     }else {
 
       G4VisAttributes* RichTbHpdEnvelopeTube_logVisAtt =
        new G4VisAttributes(G4Colour(0.2,0.2,0.8));
       if(RichTbHpdEnvelopeTube_vis == 2 ){
         RichTbHpdEnvelopeTube_logVisAtt->SetForceWireframe(true);
       }else{
         RichTbHpdEnvelopeTube_logVisAtt->SetForceWireframe(false);
         RichTbHpdEnvelopeTube_logVisAtt->SetForceSolid(true);
       }
       RichTbHpdEnvelopeTube_logVisAtt->SetVisibility(true);

       RichTbHpdEnvelopeTube_LV->
            SetVisAttributes(RichTbHpdEnvelopeTube_logVisAtt);

     }

     //hpd quartz window
      G4LogicalVolume* RichTbHpdQuartzW_LV = 
       RTbPhotoDetector->getRichHPD(ihpd)->getHpdQuartzWLogicalVolume();


     if(RichTbHpdQuartzW_vis == 0 ){
      RichTbHpdQuartzW_LV->
           SetVisAttributes (G4VisAttributes::Invisible);
     }else {
 
       G4VisAttributes* RichTbHpdQuartzW_logVisAtt =
        new G4VisAttributes(G4Colour(0.8,0.8,0.8));
       if(RichTbHpdQuartzW_vis == 2 ){
         RichTbHpdQuartzW_logVisAtt->SetForceWireframe(true);
       }else{
         RichTbHpdQuartzW_logVisAtt->SetForceWireframe(false);
         RichTbHpdQuartzW_logVisAtt->SetForceSolid(true);
       }
       RichTbHpdQuartzW_logVisAtt->SetVisibility(true);

       RichTbHpdQuartzW_LV->
            SetVisAttributes(RichTbHpdQuartzW_logVisAtt);

     }

     //hpd Ph Cathode 
      G4LogicalVolume* RichTbHpdPhCathode_LV = 
       RTbPhotoDetector->getRichHPD(ihpd)->getHpdPhCathodeLogicalVolume();


     if(RichTbHpdPhCathode_vis == 0 ){
      RichTbHpdPhCathode_LV->
           SetVisAttributes (G4VisAttributes::Invisible);
     }else {
 
       G4VisAttributes* RichTbHpdPhCathode_logVisAtt =
        new G4VisAttributes(G4Colour(0.2,0.3,0.3));
       if(RichTbHpdPhCathode_vis == 2 ){
         RichTbHpdPhCathode_logVisAtt->SetForceWireframe(true);
       }else{
         RichTbHpdPhCathode_logVisAtt->SetForceWireframe(false);
         RichTbHpdPhCathode_logVisAtt->SetForceSolid(true);
       }
       RichTbHpdPhCathode_logVisAtt->SetVisibility(true);

       RichTbHpdPhCathode_LV->
            SetVisAttributes(RichTbHpdPhCathode_logVisAtt);

     }


     //Hpd Silicon Detector sectors
     for(G4int isect=0; isect<NumberOfSiDetSectors; isect++) {
      G4LogicalVolume* RichTbHpdSiDet_LV = 
      RTbPhotoDetector->getRichHPD(ihpd)->getRichHpdSiDetSect(isect)->
               getHpdSiSectLogicalVolume();

     if(RichTbHpdSiDet_vis == 0 ){
      RichTbHpdSiDet_LV->
           SetVisAttributes (G4VisAttributes::Invisible);
     }else {
 
       G4VisAttributes* RichTbHpdSiDet_logVisAtt =
        new G4VisAttributes(G4Colour(0.6,0.6,0.9));
       if(RichTbHpdSiDet_vis == 2 ){
         RichTbHpdSiDet_logVisAtt->SetForceWireframe(true);
       }else{
         RichTbHpdSiDet_logVisAtt->SetForceWireframe(false);
         RichTbHpdSiDet_logVisAtt->SetForceSolid(true);
       }
       RichTbHpdSiDet_logVisAtt->SetVisibility(true);

       RichTbHpdSiDet_LV->
            SetVisAttributes(RichTbHpdSiDet_logVisAtt);

     }

     //Coating on the silicon sector


       G4LogicalVolume* RichTbHpdSectCoat_LV = 
       RTbPhotoDetector->getRichHPD(ihpd)->getRichHpdSiDetSect(isect)->
               getHpdSectCoatLogicalVolume();

     if(RichTbHpdSectCoat_vis == 0 ){
      RichTbHpdSectCoat_LV->
           SetVisAttributes (G4VisAttributes::Invisible);
     }else {
 
       G4VisAttributes* RichTbHpdSectCoat_logVisAtt =
        new G4VisAttributes(G4Colour(0.8,0.3,0.9));
       if(RichTbHpdSectCoat_vis == 2 ){
         RichTbHpdSectCoat_logVisAtt->SetForceWireframe(true);
       }else{
         RichTbHpdSectCoat_logVisAtt->SetForceWireframe(false);
         RichTbHpdSectCoat_logVisAtt->SetForceSolid(true);
       }
       RichTbHpdSectCoat_logVisAtt->SetVisibility(true);

       RichTbHpdSectCoat_LV->
            SetVisAttributes(RichTbHpdSectCoat_logVisAtt);

     }
      G4int numPix=RTbROGeom ->  getROPhotDet()-> getRichHPD(ihpd) ->
                                   getRichHpdSiDetSect(isect) ->getNumSiPixels();


      for(G4int ipixel=0; ipixel<numPix ; ipixel++){


	  G4LogicalVolume* RichTbHpdSiPx_LV= 
	         RTbROGeom ->  getROPhotDet() -> getRichHPD(ihpd) ->
	        getRichHpdSiDetSect(isect) -> getHpdSiPx(ipixel) ->
	        getHpdSiPixelLogicalVolume();


          if(RichTbHpdSiPx_vis == 0 ){
	    RichTbHpdSiPx_LV->SetVisAttributes (G4VisAttributes::Invisible);
	  }else {
 
	    G4VisAttributes* RichTbHpdSiPx_logVisAtt =
		    new G4VisAttributes(G4Colour(0.0,0.7,0.6));
	    if(RichTbHpdSiPx_vis == 2 ){
	      RichTbHpdSiPx_logVisAtt->SetForceWireframe(true);
	    }else{
	      RichTbHpdSiPx_logVisAtt->SetForceWireframe(false);
	      RichTbHpdSiPx_logVisAtt->SetForceSolid(true);
	    }
	    RichTbHpdSiPx_logVisAtt->SetVisibility(true);

	    RichTbHpdSiPx_LV->  SetVisAttributes(RichTbHpdSiPx_logVisAtt);

	  }




      }

     }
  }


}
RichTbGraphics::~RichTbGraphics(){;}







