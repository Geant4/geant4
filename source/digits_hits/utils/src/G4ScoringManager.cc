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
// $Id$
//

#include "G4ScoringManager.hh"
#include "G4ScoringMessenger.hh"
#include "G4ScoreQuantityMessenger.hh"
#include "G4VScoringMesh.hh"
#include "G4THitsMap.hh"
#include "G4VScoreColorMap.hh"
#include "G4DefaultLinearColorMap.hh"
#include "G4ScoreLogColorMap.hh"

G4ScoringManager* G4ScoringManager::fSManager = 0;

G4int G4ScoringManager::replicaLevel = 3;

G4ScoringManager* G4ScoringManager::GetScoringManager()
{
  if(!fSManager)
  {
    fSManager = new G4ScoringManager;
  }
  return fSManager;
}

G4ScoringManager* G4ScoringManager::GetScoringManagerIfExist()
{ return fSManager; }

G4ScoringManager::G4ScoringManager()
  : verboseLevel(0),fCurrentMesh(0)
{
  fMessenger = new G4ScoringMessenger(this);
  fQuantityMessenger = new G4ScoreQuantityMessenger(this);
  fColorMapDict = new ColorMapDict();
  fDefaultLinearColorMap = new G4DefaultLinearColorMap("defaultLinearColorMap");
  (*fColorMapDict)[fDefaultLinearColorMap->GetName()] = fDefaultLinearColorMap;
  G4VScoreColorMap * logColorMap = new G4ScoreLogColorMap("logColorMap");
  (*fColorMapDict)[logColorMap->GetName()] = logColorMap;
  writer = new G4VScoreWriter();
}

G4ScoringManager::~G4ScoringManager()
{
  if (writer) { delete writer; }
  delete fDefaultLinearColorMap;
  delete fColorMapDict;
  delete fQuantityMessenger;
  delete fMessenger;
  delete fSManager;
}

void G4ScoringManager::SetReplicaLevel(G4int lvl)
{ replicaLevel = lvl; }
G4int G4ScoringManager::GetReplicaLevel()
{ return replicaLevel; }

void G4ScoringManager::Accumulate(G4VHitsCollection* map)
{
  G4String wName = map->GetSDname();
  G4VScoringMesh* sm = FindMesh(wName);
  if(sm == 0) return;
  if(verboseLevel>9)
  { G4cout << "G4ScoringManager::Accumulate() for " << map->GetSDname() << " / " << map->GetName() << G4endl;
    G4cout << "  is calling G4VScoringMesh::Accumulate() of " << sm->GetWorldName() << G4endl; }
  sm->Accumulate(static_cast<G4THitsMap<double>*>(map));
}

G4VScoringMesh* G4ScoringManager::FindMesh(const G4String& wName)
{
  G4VScoringMesh* sm = 0;
  for(MeshVecItr itr = fMeshVec.begin(); itr != fMeshVec.end(); itr++) {
    G4String smName = (*itr)->GetWorldName();
    if(wName == smName) {
      sm = *itr;
      break;
    }
  }
  if(!sm && verboseLevel>9)
  { G4cout << "WARNING : G4ScoringManager::FindMesh() --- <" << wName << "> is not found. Null returned." << G4endl; }

  return sm;
}

void G4ScoringManager::List() const
{
  G4cout << "G4ScoringManager has " << GetNumberOfMesh() << " scoring meshes." << G4endl;
  for(MeshVecConstItr itr = fMeshVec.begin(); itr != fMeshVec.end(); itr++) {
   (*itr)->List();
  }
}

void G4ScoringManager::Dump() const
{
  for(MeshVecConstItr itr = fMeshVec.begin(); itr != fMeshVec.end(); itr++) {
   (*itr)->Dump();
  }
}

void G4ScoringManager::DrawMesh(const G4String& meshName,
                                const G4String& psName,
                                const G4String& colorMapName, G4int axflg)
{
  G4VScoringMesh* mesh = FindMesh(meshName);
  if(mesh) 
  {
    G4VScoreColorMap* colorMap = GetScoreColorMap(colorMapName);
    if(!colorMap)
    {
      G4cerr << "WARNING : Score color map <" << colorMapName << "> is not found. Default linear color map is used." << G4endl;
      colorMap = fDefaultLinearColorMap;
    }
    mesh->DrawMesh(psName,colorMap,axflg);
  } else {
    G4cerr << "ERROR : G4ScoringManager::DrawMesh() --- <"
	   << meshName << "> is not found. Nothing is done." << G4endl;
  }
}

void G4ScoringManager::DrawMesh(const G4String& meshName,
                                const G4String& psName,
                                      G4int idxPlane, G4int iColumn,
                                const G4String& colorMapName)
{
  G4VScoringMesh* mesh = FindMesh(meshName);
  if(mesh) 
  {
    G4VScoreColorMap* colorMap = GetScoreColorMap(colorMapName);
    if(!colorMap)
    {
      G4cerr << "WARNING : Score color map <" << colorMapName << "> is not found. Default linear color map is used." << G4endl;
      colorMap = fDefaultLinearColorMap;
    }
    mesh->DrawMesh(psName,idxPlane,iColumn,colorMap);
  } else {
    G4cerr << "ERROR : G4ScoringManager::DrawMesh() --- <"
	   << meshName << "> is not found. Nothing is done." << G4endl;
  }
}

void G4ScoringManager::DumpQuantityToFile(const G4String& meshName,
                                          const G4String& psName,
                                          const G4String& fileName,
                                          const G4String& option)
{
  G4VScoringMesh* mesh = FindMesh(meshName);
  if(mesh) {
    writer->SetScoringMesh(mesh);
    writer->DumpQuantityToFile(psName, fileName, option);
  } else {
    G4cerr << "ERROR : G4ScoringManager::DrawQuantityToFile() --- <"
	   << meshName << "> is not found. Nothing is done." << G4endl;
  }
}

void G4ScoringManager::DumpAllQuantitiesToFile(const G4String& meshName,
                                               const G4String& fileName,
                                               const G4String& option)
{
  G4VScoringMesh* mesh = FindMesh(meshName);
  if(mesh) {
    writer->SetScoringMesh(mesh);
    writer->DumpAllQuantitiesToFile(fileName, option);
  } else {
    G4cerr << "ERROR : G4ScoringManager::DrawAllQuantitiesToFile() --- <"
	   << meshName << "> is not found. Nothing is done." << G4endl;
  }
}

void G4ScoringManager::RegisterScoreColorMap(G4VScoreColorMap* colorMap)
{
  if(fColorMapDict->find(colorMap->GetName()) != fColorMapDict->end())
  {
    G4cerr << "ERROR : G4ScoringManager::RegisterScoreColorMap -- "
           << colorMap->GetName() << " has already been registered. Method ignored." << G4endl;
  }
  else
  {
    (*fColorMapDict)[colorMap->GetName()] = colorMap;
  }
}

G4VScoreColorMap* G4ScoringManager::GetScoreColorMap(const G4String& mapName)
{
  ColorMapDictItr mItr = fColorMapDict->find(mapName);
  if(mItr == fColorMapDict->end()) { return 0; }
  return (mItr->second);
}

void G4ScoringManager::ListScoreColorMaps()
{
  G4cout << "Registered Score Color Maps -------------------------------------------------------" << G4endl;
  ColorMapDictItr mItr = fColorMapDict->begin();
  for(;mItr!=fColorMapDict->end();mItr++)
  { G4cout << "   " << mItr->first; }
  G4cout << G4endl;
}


