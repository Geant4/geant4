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
/// \file RE06/src/RE06ParallelWorld.cc
/// \brief Implementation of the RE06ParallelWorld class
//
// 

#include "RE06ParallelWorld.hh"

#include "G4Tubs.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4PVReplica.hh"

#include "G4SDManager.hh"
#include "G4MultiFunctionalDetector.hh"
#include "G4VPrimitiveScorer.hh"
#include "G4PSEnergyDeposit.hh"
#include "G4PSNofSecondary.hh"
#include "G4PSTrackLength.hh"
#include "G4PSNofStep.hh"
#include "G4PSMinKinEAtGeneration.hh"
#include "G4VSDFilter.hh"
#include "G4SDParticleFilter.hh"
#include "G4ios.hh"
#include "G4SystemOfUnits.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4ThreadLocal G4bool RE06ParallelWorld::fSDConstructed = false;

RE06ParallelWorld::RE06ParallelWorld(G4String worldName)
:G4VUserParallelWorld(worldName),
 fConstructed(false),
 fSerial(false),
 fTotalThickness(2.0*m),
 fNumberOfLayers(20)
{
  for(size_t i=0;i<3;i++)
  {
    fCalorLogical[i] = 0;
    fLayerLogical[i] = 0;
    fCalorPhysical[i] = 0;
    fLayerPhysical[i] = 0;
  }
  fCalName[0] = "Calor-AP";
  fCalName[1] = "Calor-BP";
  fCalName[2] = "Calor-CP";
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

RE06ParallelWorld::~RE06ParallelWorld()
{;}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void RE06ParallelWorld::Construct()
{
  if(!fConstructed)
  { 
    fConstructed = true;
    SetupGeometry();
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void RE06ParallelWorld::ConstructSD()
{
  if(!fSDConstructed)
  {
    fSDConstructed = true;
    SetupDetectors();
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void RE06ParallelWorld::SetupGeometry()
{
  //     
  // World
  //
  G4VPhysicalVolume* ghostWorld = GetWorld();
  G4LogicalVolume* worldLogical = ghostWorld->GetLogicalVolume();
  
  //                               
  // Calorimeter
  //  
  G4VSolid* calorSolid 
    = new G4Tubs("Calor",0.0,0.5*m,fTotalThickness/2.,0.0,360.*deg);
  G4int i;
  for(i=0;i<3;i++)
  {
    fCalorLogical[i] = new G4LogicalVolume(calorSolid,0,fCalName[i]);
    if(fSerial)
    {
      fCalorPhysical[i] = new G4PVPlacement(0,
                 G4ThreeVector(0.,0.,G4double(i-1)*fTotalThickness),
                 fCalorLogical[i],fCalName[i],worldLogical,false,i);
    }
    else
    {
      fCalorPhysical[i] = new G4PVPlacement(0,
                 G4ThreeVector(0.,G4double(i-1)*m,0.),
                 fCalorLogical[i],fCalName[i],worldLogical,false,i);
    }
  }
 
  //                                 
  // Layers --- as absorbers
  //
  G4VSolid* layerSolid 
    = new G4Tubs("Layer",0.0,0.5*m,fTotalThickness/2.,0.0,360.*deg);
  for(i=0;i<3;i++)
  {
    fLayerLogical[i] 
      = new G4LogicalVolume(layerSolid,0,fCalName[i]+"_LayerLog");
    fLayerPhysical[i] 
      = new G4PVReplica(fCalName[i]+"_Layer",fLayerLogical[i],fCalorLogical[i],
                        kRho,fNumberOfLayers,0.5*m/fNumberOfLayers);
  }
   
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void RE06ParallelWorld::SetupDetectors()
{
  G4SDManager::GetSDMpointer()->SetVerboseLevel(1);
  G4String filterName, particleName;

  G4SDParticleFilter* gammaFilter 
    = new G4SDParticleFilter(filterName="gammaFilter",particleName="gamma");
  G4SDParticleFilter* electronFilter 
    = new G4SDParticleFilter(filterName="electronFilter",particleName="e-");
  G4SDParticleFilter* positronFilter 
    = new G4SDParticleFilter(filterName="positronFilter",particleName="e+");
  G4SDParticleFilter* epFilter 
    = new G4SDParticleFilter(filterName="epFilter");
  epFilter->add(particleName="e-");
  epFilter->add(particleName="e+");

  for(G4int i=0;i<3;i++)
  {
    G4String detName = fCalName[i]+"_para";
    G4MultiFunctionalDetector* det = new G4MultiFunctionalDetector(detName);

    G4VPrimitiveScorer* primitive;
    primitive = new G4PSEnergyDeposit("eDep");
    det->RegisterPrimitive(primitive);
    primitive = new G4PSNofSecondary("nGamma");
    primitive->SetFilter(gammaFilter);
    det->RegisterPrimitive(primitive);
    primitive = new G4PSNofSecondary("nElectron");
    primitive->SetFilter(electronFilter);
    det->RegisterPrimitive(primitive);
    primitive = new G4PSNofSecondary("nPositron");
    primitive->SetFilter(positronFilter);
    det->RegisterPrimitive(primitive);
    primitive = new G4PSTrackLength("trackLength");
    primitive->SetFilter(epFilter);
    det->RegisterPrimitive(primitive);
    primitive = new G4PSNofStep("nStep");
    primitive->SetFilter(epFilter);
    det->RegisterPrimitive(primitive);

    G4SDManager::GetSDMpointer()->AddNewDetector(det);
    SetSensitiveDetector(fLayerLogical[i], det);
  }
  G4SDManager::GetSDMpointer()->SetVerboseLevel(0);
  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void RE06ParallelWorld::SetSerialGeometry(G4bool serial)
{
  if(fSerial==serial) return;
  fSerial=serial;
  if(!fConstructed) return;
  for(G4int i=0;i<3;i++)
  {
    if(fSerial)
    { 
      fCalorPhysical[i]
        ->SetTranslation(G4ThreeVector(0.,0.,G4double(i-1)*2.*m)); 
    }
    else
    { 
      fCalorPhysical[i]
        ->SetTranslation(G4ThreeVector(0.,G4double(i-1)*m,0.)); 
    }
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
