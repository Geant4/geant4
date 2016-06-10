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
// $Id: G4tgbDetectorConstruction.cc 66363 2012-12-18 09:12:54Z gcosmo $
//
//
// class G4tgbDetectorConstruction

// History:
// - Created.                                 P.Arce, CIEMAT (November 2007)
// -------------------------------------------------------------------------

#include "G4tgbDetectorConstruction.hh"
#include "G4tgbVolume.hh"
#include "G4tgbVolumeMgr.hh"

#include "G4tgrVolume.hh"
#include "G4tgrVolumeMgr.hh"
#include "G4tgrMessenger.hh"

#include "G4Material.hh"
#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4LogicalVolume.hh"
#include "G4ThreeVector.hh"
#include "G4PVPlacement.hh"

//---------------------------------------------------------------------
G4tgbDetectorConstruction::G4tgbDetectorConstruction()
{;}

//---------------------------------------------------------------------
G4tgbDetectorConstruction::~G4tgbDetectorConstruction()
{;}

//---------------------------------------------------------------------
G4VPhysicalVolume* G4tgbDetectorConstruction::Construct()
{
  //------------------- construct g4 geometry
  //---------- find top G4tgrVolume 
  G4tgrVolumeMgr* tgrVolmgr = G4tgrVolumeMgr::GetInstance();
  const G4tgrVolume* tgrVoltop = tgrVolmgr->GetTopVolume();  

  //---------- copy list of G4tgrVolume's to list of G4tgbVolume's
  //           (just a trick to make all GEANT4 volume building in this class)
  G4tgbVolumeMgr* tgbVolmgr = G4tgbVolumeMgr::GetInstance();
  tgbVolmgr->CopyVolumes();
  //---------- find corresponding volume in list of G4tgbVolume's
  G4tgbVolume* tgbVoltop = tgbVolmgr->FindVolume( tgrVoltop->GetName() );
  
  //---------- ConstructG4Volumes of top G4tgbVolume
  //           (it will recursively build the whole tree)
  tgbVoltop->ConstructG4Volumes( 0, (const G4LogicalVolume*)0 );
 
 
  G4VPhysicalVolume* physvol = (G4tgbVolumeMgr::GetInstance())->GetTopPhysVol();

#ifdef G4VERBOSE
  if( G4tgrMessenger::GetVerboseLevel() >= 1 )
  {
    G4cout << " G4tgbDetectorConstruction::Construct() - Volume: "
           << physvol->GetName() << G4endl;
  }
#endif

  return physvol;
}
