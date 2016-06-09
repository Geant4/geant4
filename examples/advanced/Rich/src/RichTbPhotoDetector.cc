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
// RichTbPhotoDetector.cc for Rich of LHCb
// History:
// Created: Sajan Easo (Sajan.Easo@cern.ch)
// Revision and changes: Patricia Mendez (Patricia.Mendez@cern.ch)
/////////////////////////////////////////////////////////////////////////////
#include <iostream.h>
#include <fstream.h>
#include "globals.hh"
#include "RichTbPhotoDetector.hh"
#include "G4VPhysicalVolume.hh"
#include "RichTbGeometryParameters.hh"


#include "RichTbSiDet.hh"
#include "RichTbSiPixel.hh"
#include "G4ThreeVector.hh"
#include "G4RotationMatrix.hh"

RichTbPhotoDetector::RichTbPhotoDetector() {;}
RichTbPhotoDetector::RichTbPhotoDetector(RichTbMaterial* RMaterial,
  RichTbComponent* RTbComponent, RichTbRunConfig* rConfig, 
  G4bool ConstructTrackingSwitch){ 

    ConstructTrackingGeometrySwitch=ConstructTrackingSwitch;

  G4VPhysicalVolume* EnclosurePhys= 
    RTbComponent-> getEnclosurePhysicalVolume();
  //Now for the Individual HPDS

  NumOfHpds=NumberOfHpds;
  //Normally all hpds have the same number of Si sectors.
  G4int NumSect= NumberOfSiDetSectors;
  // Loop through the HPDs
  for(G4int ihpd=0; ihpd<NumberOfHpds; ihpd++){
    // create the HPDs    
    richHPD[ihpd] = new RichTbHpd(RMaterial,EnclosurePhys,ihpd,NumSect,
                                     ConstructTrackingGeometrySwitch);

  }
}
RichTbPhotoDetector::~RichTbPhotoDetector(){ ; } 










