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
//
//
// class G4tgbDetectorBuilder

// History:
// - Created.                                 P.Arce, CIEMAT (November 2007)
// -------------------------------------------------------------------------

#include "G4tgbVolume.hh"
#include "G4tgbDetectorBuilder.hh"
#include "G4tgbMaterialMgr.hh"
#include "G4tgbRotationMatrixMgr.hh"

#include "G4tgrVolumeMgr.hh"
#include "G4tgrFileReader.hh"
#include "G4tgrUtils.hh"

#include "G4VSolid.hh"
#include "G4PhysicalVolumeStore.hh"
#include "G4ReflectionFactory.hh"
#include "G4tgrMessenger.hh"
#include "G4tgbVolumeMgr.hh"


//---------------------------------------------------------------------
G4tgbDetectorBuilder::G4tgbDetectorBuilder() 
{
}

//---------------------------------------------------------------------
G4tgbDetectorBuilder::~G4tgbDetectorBuilder() 
{
}

//---------------------------------------------------------------------
const G4tgrVolume* G4tgbDetectorBuilder::ReadDetector()
{
  //------------------- construct g4 geometry
  //---------- find top G4tgrVolume 
  G4tgrFileReader::GetInstance()->ReadFiles();
  G4tgrVolumeMgr* tgrVolmgr = G4tgrVolumeMgr::GetInstance();
  const G4tgrVolume* tgrVoltop = tgrVolmgr->GetTopVolume();  
  return tgrVoltop;
}


//---------------------------------------------------------------------
G4VPhysicalVolume*
G4tgbDetectorBuilder::ConstructDetector( const G4tgrVolume* tgrVoltop )
{
  //---------- copy list of G4tgrVolume's to list of G4tgbVolume's
  //           (just a trick to make all GEANT4 volume building in this class)
  G4tgbVolumeMgr* tgbVolmgr = G4tgbVolumeMgr::GetInstance();
  tgbVolmgr->CopyVolumes();
  //---------- find corresponding volume in list of G4tgbVolume's
  G4tgbVolume* tgbVoltop = tgbVolmgr->FindVolume( tgrVoltop->GetName() );
  
  //---------- ConstructG4Volumes of top G4tgbVolume
  //           (it will recursively build the whole tree)
  tgbVoltop->ConstructG4Volumes( 0, (const G4LogicalVolume*)0 );

  G4VPhysicalVolume* physvol = tgbVolmgr->GetTopPhysVol();

#ifdef G4VERBOSE
  if( G4tgrMessenger::GetVerboseLevel() >= 1 )
  {
    G4cout << " G4tgbDetectorConstruction::ConstructDetector() - Volume: "
           << physvol->GetName() << G4endl;
  }
#endif

  return physvol;

}
