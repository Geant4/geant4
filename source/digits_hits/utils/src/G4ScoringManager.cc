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

#include "G4ScoringManager.hh"
#include "G4ScoringMessenger.hh"
#include "G4ScoreQuantityMessenger.hh"
#include "G4VScoringMesh.hh"
#include "G4THitsMap.hh"
#include "G4VScoreColorMap.hh"
#include "G4DefaultLinearColorMap.hh"
#include "G4ScoreLogColorMap.hh"

G4ThreadLocal G4ScoringManager* G4ScoringManager::fSManager = nullptr;

G4ThreadLocal G4int G4ScoringManager::replicaLevel = 3;

G4ScoringManager* G4ScoringManager::GetScoringManager()
{
  if(fSManager == nullptr)
  {
    fSManager = new G4ScoringManager;
  }
  return fSManager;
}

G4ScoringManager* G4ScoringManager::GetScoringManagerIfExist()
{
  return fSManager;
}

G4ScoringManager::G4ScoringManager()
  : verboseLevel(0)
  , fCurrentMesh(nullptr)
{
  fMessenger             = new G4ScoringMessenger(this);
  fQuantityMessenger     = new G4ScoreQuantityMessenger(this);
  fColorMapDict          = new ColorMapDict();
  fDefaultLinearColorMap = new G4DefaultLinearColorMap("defaultLinearColorMap");
  (*fColorMapDict)[fDefaultLinearColorMap->GetName()] = fDefaultLinearColorMap;
  G4VScoreColorMap* logColorMap = new G4ScoreLogColorMap("logColorMap");
  (*fColorMapDict)[logColorMap->GetName()] = logColorMap;
  writer                                   = new G4VScoreWriter();
}

G4ScoringManager::~G4ScoringManager()
{ 
  delete writer;
  delete fDefaultLinearColorMap;
  delete fColorMapDict;
  delete fQuantityMessenger;
  delete fMessenger;
  fSManager = nullptr;
}

void G4ScoringManager::SetReplicaLevel(G4int lvl) { replicaLevel = lvl; }
G4int G4ScoringManager::GetReplicaLevel() { return replicaLevel; }

void G4ScoringManager::Accumulate(G4VHitsCollection* map)
{
  auto sm = FindMesh(map);
  if(sm == nullptr)
    return;
  if(verboseLevel > 9)
  {
    G4cout << "G4ScoringManager::Accumulate() for " << map->GetSDname() << " / "
           << map->GetName() << G4endl;
    G4cout << "  is calling G4VScoringMesh::Accumulate() of "
           << sm->GetWorldName() << G4endl;
  }
  sm->Accumulate(static_cast<G4THitsMap<double>*>(map));
}

G4VScoringMesh* G4ScoringManager::FindMesh(G4VHitsCollection* map)
{
  auto colID         = map->GetColID();
  G4VScoringMesh* sm = nullptr;
  auto msh           = fMeshMap.find(colID);
  if(msh == fMeshMap.end())
  {
    auto wName      = map->GetSDname();
    sm              = FindMesh(wName);
    fMeshMap[colID] = sm;
  }
  else
  {
    sm = (*msh).second;
  }
  return sm;
}

G4VScoringMesh* G4ScoringManager::FindMesh(const G4String& wName)
{
  G4VScoringMesh* sm = nullptr;
  for(auto msh : fMeshVec)
  {
    if(msh->GetWorldName() == wName)
      return msh;
  }
  if((sm == nullptr) && verboseLevel > 9)
  {
    G4cout << "WARNING : G4ScoringManager::FindMesh() --- <" << wName
           << "> is not found. Null returned." << G4endl;
  }
  return nullptr;
}

void G4ScoringManager::List() const
{
  G4cout << "G4ScoringManager has " << GetNumberOfMesh() << " scoring meshes."
         << G4endl;
  for(auto msh : fMeshVec)
    msh->List();
}

void G4ScoringManager::Dump() const
{
  for(auto msh : fMeshVec)
    msh->Dump();
}

void G4ScoringManager::DrawMesh(const G4String& meshName,
                                const G4String& psName,
                                const G4String& colorMapName, G4int axflg)
{
  G4VScoringMesh* mesh = FindMesh(meshName);
  if(mesh != nullptr)
  {
    G4VScoreColorMap* colorMap = GetScoreColorMap(colorMapName);
    if(colorMap == nullptr)
    {
      G4cerr << "WARNING : Score color map <" << colorMapName
             << "> is not found. Default linear color map is used." << G4endl;
      colorMap = fDefaultLinearColorMap;
    }
    mesh->DrawMesh(psName, colorMap, axflg);
  }
  else
  {
    G4cerr << "ERROR : G4ScoringManager::DrawMesh() --- <" << meshName
           << "> is not found. Nothing is done." << G4endl;
  }
}

void G4ScoringManager::DrawMesh(const G4String& meshName,
                                const G4String& psName, G4int idxPlane,
                                G4int iColumn, const G4String& colorMapName)
{
  G4VScoringMesh* mesh = FindMesh(meshName);
  if(mesh != nullptr)
  {
    G4VScoreColorMap* colorMap = GetScoreColorMap(colorMapName);
    if(colorMap == nullptr)
    {
      G4cerr << "WARNING : Score color map <" << colorMapName
             << "> is not found. Default linear color map is used." << G4endl;
      colorMap = fDefaultLinearColorMap;
    }
    mesh->DrawMesh(psName, idxPlane, iColumn, colorMap);
  }
  else
  {
    G4cerr << "ERROR : G4ScoringManager::DrawMesh() --- <" << meshName
           << "> is not found. Nothing is done." << G4endl;
  }
}

void G4ScoringManager::DumpQuantityToFile(const G4String& meshName,
                                          const G4String& psName,
                                          const G4String& fileName,
                                          const G4String& option)
{
  G4VScoringMesh* mesh = FindMesh(meshName);
  if(mesh != nullptr)
  {
    writer->SetScoringMesh(mesh);
    writer->DumpQuantityToFile(psName, fileName, option);
  }
  else
  {
    G4cerr << "ERROR : G4ScoringManager::DrawQuantityToFile() --- <" << meshName
           << "> is not found. Nothing is done." << G4endl;
  }
}

void G4ScoringManager::DumpAllQuantitiesToFile(const G4String& meshName,
                                               const G4String& fileName,
                                               const G4String& option)
{
  G4VScoringMesh* mesh = FindMesh(meshName);
  if(mesh != nullptr)
  {
    writer->SetScoringMesh(mesh);
    writer->DumpAllQuantitiesToFile(fileName, option);
  }
  else
  {
    G4cerr << "ERROR : G4ScoringManager::DrawAllQuantitiesToFile() --- <"
           << meshName << "> is not found. Nothing is done." << G4endl;
  }
}

void G4ScoringManager::RegisterScoreColorMap(G4VScoreColorMap* colorMap)
{
  if(fColorMapDict->find(colorMap->GetName()) != fColorMapDict->cend())
  {
    G4cerr << "ERROR : G4ScoringManager::RegisterScoreColorMap -- "
           << colorMap->GetName()
           << " has already been registered. Method ignored." << G4endl;
  }
  else
  {
    (*fColorMapDict)[colorMap->GetName()] = colorMap;
  }
}

G4VScoreColorMap* G4ScoringManager::GetScoreColorMap(const G4String& mapName)
{
  auto mItr = fColorMapDict->find(mapName);
  if(mItr == fColorMapDict->cend())
  {
    return nullptr;
  }
  return (mItr->second);
}

void G4ScoringManager::ListScoreColorMaps()
{
  G4cout << "Registered Score Color Maps "
            "-------------------------------------------------------"
         << G4endl;
  auto mItr = fColorMapDict->cbegin();
  for(; mItr != fColorMapDict->cend(); ++mItr)
  {
    G4cout << "   " << mItr->first;
  }
  G4cout << G4endl;
}

void G4ScoringManager::Merge(const G4ScoringManager* mgr)
{
  for(G4int i = 0; i < (G4int)GetNumberOfMesh(); ++i)
  {
    G4VScoringMesh* fMesh  = GetMesh(i);
    G4VScoringMesh* scMesh = mgr->GetMesh(i);
    fMesh->Merge(scMesh);
  }
}
