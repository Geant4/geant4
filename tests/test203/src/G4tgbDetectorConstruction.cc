
#include "G4tgbDetectorConstruction.hh"
#include "G4tgbVolume.hh"
#include "G4tgbVolumeMgr.hh"

#include "G4tgrVolume.hh"
#include "G4tgrVolumeMgr.hh"

#include "G4Material.hh"
#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4LogicalVolume.hh"
#include "G4ThreeVector.hh"
#include "G4PVPlacement.hh"
#include "globals.hh"

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

  //---------- copy list of G4tgrVolume's to list of G4tgbVolume's (just a trick to make all GEANT4 volume building in this class)
  G4tgbVolumeMgr* tgbVolmgr = G4tgbVolumeMgr::GetInstance();
  tgbVolmgr->CopyVolumes();
  //---------- find corresponding volume in list of G4tgbVolume's
  G4tgbVolume* tgbVoltop = tgbVolmgr->FindVolume( tgrVoltop->GetName() );
  
  //---------- ConstructG4Volumes of top G4tgbVolume (it will recursively build the whole tree)
  tgbVoltop->ConstructG4Volumes( 0, (const G4LogicalVolume*)0 );
 
 
  G4VPhysicalVolume* physvol = (G4tgbVolumeMgr::GetInstance())->GetTopPhysVol();
  cout << " G4tgbDetectorConstruction::Construct. TopPV " << physvol->GetName() << endl;
  return physvol;
    //  return (G4tgbVolumeMgr::GetInstance())->getTopPV();
}

