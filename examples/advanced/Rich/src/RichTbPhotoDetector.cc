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
// RichTbPhotoDetector.cc for Rich of LHCb
// History:
// Created: Sajan Easo (Sajan.Easo@cern.ch)
// Revision and changes: Patricia Mendez (Patricia.Mendez@cern.ch)
/////////////////////////////////////////////////////////////////////////////
#include <iostream>
#include <fstream>
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
  RichTbComponent* RTbComponent, RichTbRunConfig*, 
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










