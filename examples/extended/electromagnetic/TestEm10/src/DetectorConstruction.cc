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
/// \file electromagnetic/TestEm10/src/DetectorConstruction.cc
/// \brief Implementation of the DetectorConstruction class
//
//
//
//

#include "DetectorConstruction.hh"
#include "DetectorSimpleALICE.hh"
#include "DetectorALICE06.hh"
#include "DetectorBari05.hh"
#include "DetectorHarris73.hh"
#include "DetectorWatase86.hh"
#include "DetectorBarr90.hh"
#include "DetectorMessenger.hh"
#include "Materials.hh"
#include "G4Region.hh"
#include "G4RegionStore.hh"
#include "G4ProductionCuts.hh"

#include "G4LogicalVolumeStore.hh"
#include "G4LogicalVolume.hh"
#include "G4ios.hh"

#include "G4SystemOfUnits.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetectorConstruction::DetectorConstruction()
  :G4VUserDetectorConstruction(),
  fRadiatorDescription(0),
  fSetUp("simpleALICE")
{
  fDetectorMessenger = new DetectorMessenger(this);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetectorConstruction::~DetectorConstruction()
{
  delete fRadiatorDescription;
  delete fDetectorMessenger;
  delete Materials::GetInstance();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume* DetectorConstruction::Construct()
{
  G4VPhysicalVolume* world = nullptr;
  if( fSetUp == "simpleALICE" )
  {
    DetectorSimpleALICE  simpleALICE;
    world = simpleALICE.Construct();
    fRadiatorDescription = simpleALICE.GetRadiatorDescription();
  }
  else if( fSetUp == "alice06" )
  {
    DetectorALICE06  alice06;
    world = alice06.Construct();
    fRadiatorDescription = alice06.GetRadiatorDescription();
  }
  else if( fSetUp == "bari05" )
  {
    DetectorBari05  bari05;
    world = bari05.Construct();
    fRadiatorDescription = bari05.GetRadiatorDescription();
  }
  else if( fSetUp == "harris73" )
  {
    DetectorHarris73  harris73;
    world = harris73.Construct();
    fRadiatorDescription = harris73.GetRadiatorDescription();
  }
  else if( fSetUp == "watase86" )
  {
    DetectorWatase86  watase86;
    world = watase86.Construct();
    fRadiatorDescription = watase86.GetRadiatorDescription();
  }
  else if( fSetUp == "barr90" )
  {
    DetectorBarr90  barr90;
    world = barr90.Construct();
    fRadiatorDescription = barr90.GetRadiatorDescription();
  }
  else
  {
    G4cout << "Experimental setup is unsupported. Check /XTRdetector/setup " 
           <<G4endl;
    G4cout << "Run default: simpleALICE"<<G4endl;

    DetectorSimpleALICE  simpleALICE;
    world = simpleALICE.Construct();
    fRadiatorDescription = simpleALICE.GetRadiatorDescription();
  }
  G4RegionStore* regionStore = G4RegionStore::GetInstance();
  size_t nRegions = regionStore->size();
  G4double cut = 1.*CLHEP::mm;
  for(size_t i=0; i<nRegions; ++i) {
    G4Region* reg = (*regionStore)[i];
    if(!reg->GetProductionCuts()) {
      G4ProductionCuts* cuts = new G4ProductionCuts();
      cuts->SetProductionCut(cut,"gamma");
      cuts->SetProductionCut(cut,"e-");
      cuts->SetProductionCut(cut,"e+");
      cuts->SetProductionCut(cut,"proton");
      reg->SetProductionCuts(cuts);
    }
  }
  return world;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

RadiatorDescription* DetectorConstruction::GetRadiatorDescription() const
{
  if ( ! fRadiatorDescription ) {
    G4cout << "RadiatorDescription is not defined" << G4endl;
  }
  return fRadiatorDescription;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4Material* DetectorConstruction::GetAbsorberMaterial() const
{
  G4LogicalVolume* lv
    = G4LogicalVolumeStore::GetInstance()->GetVolume("Absorber");

  if ( ! lv ) {
    G4cerr << "Absorber logical volume is not defined." << G4endl;
    return 0;
  }

  return lv->GetMaterial();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

