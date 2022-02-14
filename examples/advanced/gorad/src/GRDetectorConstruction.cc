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
// GRDetectorConstruction.cc
//   Read a GDML file to set up the geometry.
//   This class also takes care of several utility methods on geometry
//   and creates a parallel world for geometry importance biasing if
//   requested.
//
// History
//   September 8th, 2020 : first implementation
//
// ********************************************************************

#include "GRDetectorConstruction.hh"

#include "G4GDMLParser.hh"
#include "GRDetectorConstructionMessenger.hh"
#include "GRGeomImpBiasWorld.hh"

#include "G4VPhysicalVolume.hh"
#include "G4RunManager.hh"
#include "G4LogicalVolume.hh"
#include "G4Box.hh"
#include <algorithm>

G4double GRDetectorConstruction::worldSize = -1.;

GRDetectorConstruction::GRDetectorConstruction()
: gdmlFile("noName"), fWorld(nullptr), initialized(false)
{
  messenger = new GRDetectorConstructionMessenger(this);
  parser = new G4GDMLParser();
}

GRDetectorConstruction::~GRDetectorConstruction()
{
  delete messenger;
  delete parser;
}

G4VPhysicalVolume* GRDetectorConstruction::Construct()
{
  if(!initialized)
  {
    Read();
    if(applyGeomImpBias)
    {
      G4String wName = "GeomBias";
      geomImpBiasWorld = new GRGeomImpBiasWorld(wName,this);
      RegisterParallelWorld(geomImpBiasWorld);
    }
  }
  return fWorld;
}

void GRDetectorConstruction::ConstructSDAndField()
{ ; }

G4bool GRDetectorConstruction::SetGDMLFile(G4String& gdml)
{
  G4bool valid = true; // parser->IsValid(gdml); ## GDML parser fix needed
  if(valid)
  {
    if(initialized)
    {
      parser->Clear();
      G4RunManager::GetRunManager()->ReinitializeGeometry(true);
    }
    gdmlFile = gdml;
  }
  return valid;
}

void GRDetectorConstruction::Read()
{
  parser->Read(gdmlFile);
  fWorld = parser->GetWorldVolume();
  const G4Box* worldBox = dynamic_cast<const G4Box*>(fWorld->GetLogicalVolume()->GetSolid());
  if(!worldBox)
  {
    G4ExceptionDescription ed;
    ed << "PANIC!!! World volume defined in "<<gdmlFile<<" is not Box!!!";
    G4Exception("GRDetectorConstruction::Read()","GOoradGeom00012",FatalException,ed);
  }
  worldSize = std::min( {worldBox->GetXHalfLength(), worldBox->GetYHalfLength(), worldBox->GetZHalfLength()},
                        [](G4double a,G4double b) {return a<b;} );
  initialized = true;
}

#include "G4UnitsTable.hh"
#include "G4VSolid.hh"
#include "G4SolidStore.hh"
void GRDetectorConstruction::ListSolids(G4int lvl)
{
  G4cout << "*********** List of registered solids *************" << G4endl;
  auto store = G4SolidStore::GetInstance();
  auto itr = store->begin();
  for(;itr!=store->end();itr++)
  {
    switch(lvl)
    {
      case 0:
        G4cout << (*itr)->GetName() << G4endl;
        break;
      case 1:
        G4cout << (*itr)->GetName()
               << "\t volume = " << G4BestUnit((*itr)->GetCubicVolume(),"Volume")
               << "\t surface = " << G4BestUnit((*itr)->GetSurfaceArea(),"Surface")
               << G4endl;
        break;
      default:
        (*itr)->DumpInfo();
        break;
    }
  }
}

#include "G4LogicalVolume.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4Material.hh"
#include "G4VSensitiveDetector.hh"
void GRDetectorConstruction::ListLogVols(G4int lvl)
{
  G4cout << "*********** List of registered logical volumes *************" << G4endl;
  auto store = G4LogicalVolumeStore::GetInstance();
  auto itr = store->begin();
  for(;itr!=store->end();itr++)
  {
    G4cout << (*itr)->GetName() << "\t Solid = " << (*itr)->GetSolid()->GetName();
    if((*itr)->GetMaterial())
    { G4cout << "\t Material = " << (*itr)->GetMaterial()->GetName() << G4endl; }
    else
    { G4cout << "\t Material : not defined " << G4endl; }
    if(lvl<1) continue;
    G4cout << "\t region = ";
    if((*itr)->GetRegion())
    { G4cout << (*itr)->GetRegion()->GetName(); }
    else
    { G4cout << "not defined"; }
    G4cout << "\t sensitive detector = ";
    if((*itr)->GetSensitiveDetector())
    { G4cout << (*itr)->GetSensitiveDetector()->GetName(); }
    else
    { G4cout << "not defined"; }
    G4cout << G4endl;
    G4cout << "\t daughters = " << (*itr)->GetNoDaughters();
    if((*itr)->GetNoDaughters()>0)
    {
      switch((*itr)->CharacteriseDaughters())
      {
        case kNormal:
          G4cout << " (placement)"; break;
        case kReplica:
          G4cout << " (replica : " << (*itr)->GetDaughter(0)->GetMultiplicity() << ")"; break;
        case kParameterised:
          G4cout << " (parameterized : " << (*itr)->GetDaughter(0)->GetMultiplicity() << ")"; break;
        default:
          ;
      }
    }
    G4cout << G4endl;
    if(lvl<2) continue;
    if((*itr)->GetMaterial())
    { G4cout << "\t weight = " << G4BestUnit((*itr)->GetMass(),"Mass") << G4endl; }
    else
    { G4cout << "\t weight : not available" << G4endl; }
  }
}

#include "G4VPhysicalVolume.hh"
#include "G4PhysicalVolumeStore.hh"
void GRDetectorConstruction::ListPhysVols(G4int lvl)
{
  G4cout << "*********** List of registered physical volumes *************" << G4endl;
  auto store = G4PhysicalVolumeStore::GetInstance();
  auto itr = store->begin();
  for(;itr!=store->end();itr++)
  {
    switch(lvl)
    {
      case 0:
        G4cout << (*itr)->GetName() << G4endl;
        break;
      case 1:
        G4cout << (*itr)->GetName()
               << "\t logical volume = " << (*itr)->GetLogicalVolume()->GetName()
               << "\t mother logical = ";
        if((*itr)->GetMotherLogical())
        { G4cout << (*itr)->GetMotherLogical()->GetName(); }
        else
        { G4cout << "not defined"; }
        G4cout << G4endl;
        break;
      default:
        G4cout << (*itr)->GetName()
               << "\t logical volume = " << (*itr)->GetLogicalVolume()->GetName()
               << "\t mother logical = ";
        if((*itr)->GetMotherLogical())
        { G4cout << (*itr)->GetMotherLogical()->GetName(); }
        else
        { G4cout << "not defined"; }
        G4cout << "\t type = ";
        switch((*itr)->VolumeType())
        {
          case kNormal:
            G4cout << "placement"; break;
          case kReplica:
            G4cout << "replica"; break;
          case kParameterised:
            G4cout << "parameterized"; break;
          default:
            ;
        }
        G4cout << G4endl;
    }
  }
}

G4bool GRDetectorConstruction::CheckOverlap(G4String& physVolName, G4int nSpots,
                 G4int maxErr, G4double tol)
{
  G4cout << "*********** Checking overlap for <" << physVolName << "> *************" << G4endl;
  G4bool checkAll = (physVolName=="**ALL**");
  auto store = G4PhysicalVolumeStore::GetInstance();
  auto itr = store->begin();
  G4VPhysicalVolume* physVol = nullptr;
  for(;itr!=store->end();itr++)
  {
    if(checkAll || (*itr)->GetName()==physVolName)
    {
      physVol = (*itr);
      physVol->CheckOverlaps(nSpots,tol,true,maxErr);
      if(!checkAll) break;
    }
  }
  return (physVol!=nullptr);
}

#include "G4Region.hh"
#include "G4RegionStore.hh"
#include "G4RunManagerKernel.hh"
void GRDetectorConstruction::ListRegions(G4int lvl)
{
  if(lvl==2)
  {
    G4RunManagerKernel::GetRunManagerKernel()->DumpRegion();
    return;
  }
  G4cout << "*********** List of registered regions *************" << G4endl;
  auto store = G4RegionStore::GetInstance();
  auto itr = store->begin();
  for(;itr!=store->end();itr++)
  {
    G4cout << (*itr)->GetName();
    if((*itr)->GetWorldPhysical())
    {
      G4cout << "\t in the world volume <" << (*itr)->GetWorldPhysical()->GetName() << "> ";
      if((*itr)->IsInMassGeometry()) G4cout << "-- mass world";
      if((*itr)->IsInParallelGeometry()) G4cout << "-- parallel world";
    }
    else
    { G4cout << " -- is not associated to any world."; }
    G4cout << G4endl;
    if(lvl==0) continue;
    G4cout << "\t\t Root logical volume(s) : ";
    size_t nRootLV = (*itr)->GetNumberOfRootVolumes();
    std::vector<G4LogicalVolume*>::iterator lvItr = (*itr)->GetRootLogicalVolumeIterator();
    for(size_t j=0;j<nRootLV;j++)
    { G4cout << (*lvItr)->GetName() << " "; lvItr++; }
    G4cout << G4endl;
    G4cout << "\t\t Pointers : G4VUserRegionInformation[" << (*itr)->GetUserInformation()
           << "], G4UserLimits[" << (*itr)->GetUserLimits()
           << "], G4FastSimulationManager[" << (*itr)->GetFastSimulationManager()
           << "], G4UserSteppingAction[" << (*itr)->GetRegionalSteppingAction() << "]" << G4endl;
  }
}

G4bool GRDetectorConstruction::CreateRegion(G4String& regionName,G4String& logVolName)
{
  auto logVolStore = G4LogicalVolumeStore::GetInstance();
  auto itr = logVolStore->begin();
  G4LogicalVolume* logVol = nullptr;
  for(;itr!=logVolStore->end();itr++)
  {
    if((*itr)->GetName() == logVolName)
    { logVol = (*itr); break; }
  }
  if(!logVol) return false;

  auto regionStore = G4RegionStore::GetInstance();
  auto region = regionStore->FindOrCreateRegion(regionName);
  logVol->SetRegion(region);
  region->AddRootLogicalVolume(logVol);
  return true;
}
  
#include "G4MaterialTable.hh"
void GRDetectorConstruction::ListAllMaterial()
{
  auto materialTable = G4Material::GetMaterialTable();
  auto matItr = materialTable->begin();
  G4cout << "*********** List of instantiated materials **************" << G4endl;
  G4int i = 0;
  for(;matItr!=materialTable->end();matItr++)
  {
    G4cout << (*matItr)->GetName() << "\t";
    if(++i%5==0) G4cout << G4endl;
  }
  G4cout << G4endl;
}

G4bool GRDetectorConstruction::ListMaterial(G4String& matName)
{
  auto materialTable = G4Material::GetMaterialTable();
  auto matItr = materialTable->begin();
  for(;matItr!=materialTable->end();matItr++)
  {
    if((*matItr)->GetName()==matName) 
    {
      G4cout << *matItr << G4endl;
      return true;
    }
  }
  return false;
}

#include "G4NistManager.hh"
void GRDetectorConstruction::DumpNistMaterials()
{
  auto nameVec = G4NistManager::Instance()->GetNistMaterialNames();
  auto itr = nameVec.begin();
  G4int i = 0;
  for(;itr!=nameVec.end();itr++)
  {
    G4cout << std::setw(26) << *itr;
    if(++i%3==0) G4cout << G4endl;
  }
  G4cout << G4endl;
}

G4bool GRDetectorConstruction::CreateMaterial(G4String& matName)
{
  auto mat = G4NistManager::Instance()->FindOrBuildMaterial(matName);
  return (mat!=nullptr);
}

G4bool GRDetectorConstruction::GetMaterial(G4String& logVol)
{
  auto store = G4LogicalVolumeStore::GetInstance();
  std::vector<G4LogicalVolume*>::iterator itr = store->begin();
  for(;itr!=store->end();itr++)
  {
    if((*itr)->GetName()==logVol)
    {
      G4cout << "Logical volume <" << (*itr)->GetName() << "> is made of <"
           << (*itr)->GetMaterial()->GetName() << ">" << G4endl;
      return true;
    }
  }
  return false;
}

G4int GRDetectorConstruction::SetMaterial(G4String& logVolName,G4String& matName)
{
  G4LogicalVolume* logVol = nullptr;
  G4Material* mat = nullptr;

  auto store = G4LogicalVolumeStore::GetInstance();
  auto itr = store->begin();
  for(;itr!=store->end();itr++)
  {
    if((*itr)->GetName()==logVolName)
    {
      logVol = *itr;
      break;
    }
  }

  auto materialTable = G4Material::GetMaterialTable();
  auto matItr = materialTable->begin();
  for(;matItr!=materialTable->end();matItr++)
  {
    if((*matItr)->GetName()==matName) 
    {
      mat = *matItr;
      break;
    }
  }

  G4int retVal = 0;
  if(!logVol && !mat) 
  { retVal = 3; }
  else if(!logVol)
  { retVal = 1; }
  else if(!mat)
  { retVal = 2; }
  else
  { logVol->SetMaterial(mat); }

  return retVal;
}


