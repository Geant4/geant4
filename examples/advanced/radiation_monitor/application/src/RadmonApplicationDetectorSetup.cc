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
// File name:     RadmonApplicationDetectorSetup.cc
// Creation date: Sep 2005
// Main author:   Riccardo Capra <capra@ge.infn.it>
//
// Id:            $Id: RadmonApplicationDetectorSetup.cc,v 1.1.2.2 2006/06/29 16:08:33 gunter Exp $
// Tag:           $Name: geant4-08-01 $
//

// Include files
#include "RadmonApplicationDetectorSetup.hh"
#include "RadmonApplicationOptions.hh"

#include "RadmonDetectorLabelledEntitiesConstructorsFactory.hh"

#include "RadmonDetectorFlatVolumeConstructor.hh"
#include "RadmonDetectorFlatVolumeWithHoleConstructor.hh"
#include "RadmonDetectorFlatVolumeWithGroundConstructor.hh"
#include "RadmonDetectorFlatVolumeWithGroundAndKeyMarksConstructor.hh"
#include "RadmonDetectorFlatVolumeWithGroundAndKeyMarksAndHoleConstructor.hh"
#include "RadmonDetectorFlatVolumeWithTracksConstructor.hh"
#include "RadmonDetectorFlatVolumeWithPinsConstructor.hh"
#include "RadmonDetectorFlatVolumeWithPadsConstructor.hh"
#include "RadmonDetectorFlatVolumeWithTracksAndHoleConstructor.hh"
#include "RadmonDetectorCarvedFlatVolumeWithHoleConstructor.hh"
#include "RadmonDetectorCarvedFlatVolumeConstructor.hh"
#include "RadmonDetectorCarvedFlatVolumeWithGroundConstructor.hh"
#include "RadmonDetectorCarvedFlatVolumeWithGroundAndKeyMarksConstructor.hh"
#include "RadmonDetectorCarvedFlatVolumeWithGroundAndKeyMarksAndHoleConstructor.hh"
#include "RadmonDetectorCarvedFlatVolumeWithTracksConstructor.hh"
#include "RadmonDetectorCarvedFlatVolumeWithPinsConstructor.hh"
#include "RadmonDetectorCarvedFlatVolumeWithTracksAndHoleConstructor.hh"



#define DECLARE_DETECTOR_CONSTRUCTOR(name)      constructor=new name();                                                                  \
                                                if (constructor==0)                                                                      \
                                                {                                                                                        \
                                                 G4cerr << currentOptions.ApplicationName() << ": Cannot allocate " #name "." << G4endl; \
                                                 return false;                                                                           \
                                                }                                                                                        \
                                                factory->AppendLabelledEntityConstructor(constructor)

G4bool                                          RadmonApplicationDetectorSetup :: CreateDetectorEntityConstructors(RadmonDetectorLabelledEntitiesConstructorsFactory * factory)
{
 RadmonVDetectorLabelledEntityConstructor * constructor;
 
 DECLARE_DETECTOR_CONSTRUCTOR(RadmonDetectorFlatVolumeConstructor);
 DECLARE_DETECTOR_CONSTRUCTOR(RadmonDetectorFlatVolumeWithHoleConstructor);
 DECLARE_DETECTOR_CONSTRUCTOR(RadmonDetectorFlatVolumeWithGroundConstructor);
 DECLARE_DETECTOR_CONSTRUCTOR(RadmonDetectorFlatVolumeWithGroundAndKeyMarksConstructor);
 DECLARE_DETECTOR_CONSTRUCTOR(RadmonDetectorFlatVolumeWithGroundAndKeyMarksAndHoleConstructor);
 DECLARE_DETECTOR_CONSTRUCTOR(RadmonDetectorFlatVolumeWithTracksConstructor);
 DECLARE_DETECTOR_CONSTRUCTOR(RadmonDetectorFlatVolumeWithPinsConstructor);
 DECLARE_DETECTOR_CONSTRUCTOR(RadmonDetectorFlatVolumeWithTracksAndHoleConstructor);
 DECLARE_DETECTOR_CONSTRUCTOR(RadmonDetectorFlatVolumeWithPadsConstructor);
 DECLARE_DETECTOR_CONSTRUCTOR(RadmonDetectorCarvedFlatVolumeConstructor);
 DECLARE_DETECTOR_CONSTRUCTOR(RadmonDetectorCarvedFlatVolumeWithHoleConstructor);
 DECLARE_DETECTOR_CONSTRUCTOR(RadmonDetectorCarvedFlatVolumeWithTracksConstructor);
 DECLARE_DETECTOR_CONSTRUCTOR(RadmonDetectorCarvedFlatVolumeWithPinsConstructor);
 DECLARE_DETECTOR_CONSTRUCTOR(RadmonDetectorCarvedFlatVolumeWithTracksAndHoleConstructor);
 DECLARE_DETECTOR_CONSTRUCTOR(RadmonDetectorCarvedFlatVolumeWithGroundConstructor);
 DECLARE_DETECTOR_CONSTRUCTOR(RadmonDetectorCarvedFlatVolumeWithGroundAndKeyMarksConstructor);
 DECLARE_DETECTOR_CONSTRUCTOR(RadmonDetectorCarvedFlatVolumeWithGroundAndKeyMarksAndHoleConstructor);
 
 return true;
}
