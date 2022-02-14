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
//  Gorad (Geant4 Open-source Radiation Analysis and Design)
//
//  Author : Makoto Asai (SLAC National Accelerator Laboratory)
//
//  Development of Gorad is funded by NASA Johnson Space Center (JSC)
//  under the contract NNJ15HK11B.
//
// ********************************************************************
//
// GRGeomImpBiasWorld.cc
//   A parallel world class that defines the geometry
//   of the geometry improtance biasing.
//
// History
//   September 8th, 2020 : first implementation
//
// ********************************************************************

#include "GRGeomImpBiasWorld.hh"
#include "GRDetectorConstruction.hh"

#include "G4Orb.hh"
#include "G4LogicalVolume.hh"
#include "G4VPhysicalVolume.hh"
#include "G4PVPlacement.hh"

#include "G4Region.hh"
#include "GRBiasingRegionInfo.hh"

#include "G4UIcommand.hh"
#include "G4VisAttributes.hh"

GRGeomImpBiasWorld::GRGeomImpBiasWorld(G4String& wName,GRDetectorConstruction* det)
: G4VUserParallelWorld(wName),detector(det)
{ 
  stateNotifier = new GRGeomImpBiasWorldStateNotifier(this);
}

GRGeomImpBiasWorld::~GRGeomImpBiasWorld()
{ 
  delete stateNotifier;
}

void GRGeomImpBiasWorld::Construct()
{
  if(constructed) return;

  constructed = true;
  fWorld = GetWorld(); //physical volume of the world volume of parallel world
  G4LogicalVolume* motherLog = fWorld->GetLogicalVolume();
  biasingRegion = new G4Region(fWorldName+"_Region");
  auto biasInfo = new GRBiasingRegionInfo();
  biasingRegion->SetUserInformation(biasInfo);
  biasingRegion->AddRootLogicalVolume(motherLog);
  biasingRegion->SetWorld(fWorld);
  G4VisAttributes* wvisatt = new G4VisAttributes(G4Colour(.2,.2,0.));
  wvisatt->SetVisibility(false);
  motherLog->SetVisAttributes(wvisatt);

  G4double r0 = detector->geoImpP.radius;
  if(r0 < 0.)
  {
    // radius of the outermost sphere is not specified, so seting it to the default
    // value as the 80% of the world volume.
    r0 = GRDetectorConstruction::GetWorldSize() * 0.8;
    G4cout<<"############ Radius of the outermost biasing sphere is set to "<<r0<<" (mm)"<<G4endl;
  }

  G4int nL = detector->geoImpP.nLayer;
  G4ThreeVector dp = (detector->geoImpP.posT - detector->geoImpP.pos0) / (nL-1);
  G4double rt = detector->geoImpP.radiusT;
  if(rt<0.) 
  { rt = r0/nL; }
  G4double dr = (r0 - rt)/(nL-1);

  for(G4int i=0;i<nL;i++)
  {
    G4double r = r0 - dr*i;
    G4String vName = fWorldName + "_" + G4UIcommand::ConvertToString(i);
    auto sph = new G4Orb(vName+"_solid",r);
    auto lv = new G4LogicalVolume(sph,nullptr,vName+"_lv");
    G4ThreeVector pos = detector->geoImpP.pos0;
    if(i!=0) pos = dp;
    new G4PVPlacement(0,pos,lv,vName+"_pv",motherLog,false,i+1);
    motherLog = lv;
    G4VisAttributes* visatt = new G4VisAttributes(G4Colour(.7,.7,0.));
    visatt->SetVisibility(true);
    lv->SetVisAttributes(visatt);
  }
  
  return;
}

void GRGeomImpBiasWorld::ConstructSD()
{ ; }

void GRGeomImpBiasWorld::Update()
{
  GRBiasingRegionInfo* biasInfo = static_cast<GRBiasingRegionInfo*>(biasingRegion->GetUserInformation());
  biasInfo->SetBiasingFactor(detector->geoImpP.factor);
  biasInfo->SetProbability(detector->geoImpP.prob);
}


