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

#ifndef G4ScoringManager_h
#define G4ScoringManager_h 1

#include "globals.hh"
#include "G4VScoringMesh.hh"
#include <vector>
#include <map>
class G4ScoringMessenger;
class G4ScoreQuantityMessenger;
class G4VHitsCollection;
class G4VScoreColorMap;
#include "G4VScoreWriter.hh"

// class description:
//
//  This is a singleton class which manages the interactive scoring.
// The user cannot access to the constructor. The pointer of the
// only existing object can be got via G4ScoringManager::GetScoringManager()
// static method. The first invokation of this static method makes
// the singleton object.
//

using MeshVec = std::vector<G4VScoringMesh *>;
using MeshVecItr = MeshVec::iterator;
using MeshVecConstItr = MeshVec::const_iterator;
using ColorMapDict = std::map<G4String, G4VScoreColorMap *>;
using ColorMapDictItr = ColorMapDict::iterator;
using ColorMapDictConstItr = ColorMapDict::const_iterator;
using MeshMap = std::map<G4int, G4VScoringMesh *>;
using MeshMapItr = MeshMap::iterator;
using MeshMapConstItr = MeshMap::const_iterator;

class G4ScoringManager
{
 public:  // with description
  static G4ScoringManager* GetScoringManager();
  // Returns the pointer to the singleton object.
  
  static G4ScoringManager* GetScoringManagerIfExist();

  ~G4ScoringManager();

 public:
  static void SetReplicaLevel(G4int);
  static G4int GetReplicaLevel();

  void RegisterScoreColorMap(G4VScoreColorMap* colorMap);
  // Register a color map. Once registered, it is available by /score/draw and
  // /score/drawColumn commands.

  void Accumulate(G4VHitsCollection* map);
  void Merge(const G4ScoringManager* scMan);
  G4VScoringMesh* FindMesh(G4VHitsCollection* map);
  G4VScoringMesh* FindMesh(const G4String&);
  void List() const;
  void Dump() const;
  void DrawMesh(const G4String& meshName, const G4String& psName,
                const G4String& colorMapName, G4int axflg = 111);
  void DrawMesh(const G4String& meshName, const G4String& psName,
                G4int idxPlane, G4int iColumn, const G4String& colorMapName);
  void DumpQuantityToFile(const G4String& meshName, const G4String& psName,
                          const G4String& fileName,
                          const G4String& option = "");
  void DumpAllQuantitiesToFile(const G4String& meshName,
                               const G4String& fileName,
                               const G4String& option = "");
  G4VScoreColorMap* GetScoreColorMap(const G4String& mapName);
  void ListScoreColorMaps();

  inline void SetCurrentMesh(G4VScoringMesh* scm) { fCurrentMesh = scm; }
  inline G4VScoringMesh* GetCurrentMesh() const { return fCurrentMesh; }
  inline void CloseCurrentMesh() { fCurrentMesh = nullptr; }
  inline void SetVerboseLevel(G4int vl)
  {
    verboseLevel = vl;
    for(auto& itr : fMeshVec)
    {
      itr->SetVerboseLevel(vl);
    }
    if(writer != nullptr)
      writer->SetVerboseLevel(vl);
  }
  inline G4int GetVerboseLevel() const { return verboseLevel; }
  inline size_t GetNumberOfMesh() const { return fMeshVec.size(); }
  inline void RegisterScoringMesh(G4VScoringMesh* scm)
  {
    scm->SetVerboseLevel(verboseLevel);
    fMeshVec.push_back(scm);
    SetCurrentMesh(scm);
  }
  inline G4VScoringMesh* GetMesh(G4int i) const { return fMeshVec[i]; }
  inline G4String GetWorldName(G4int i) const
  {
    return fMeshVec[i]->GetWorldName();
  }

 public:  // with description
  inline void SetScoreWriter(G4VScoreWriter* sw)
  {
    delete writer;
    writer = sw;
    if(writer != nullptr)
      writer->SetVerboseLevel(verboseLevel);
  }
  // Replace score writers.

 public:
  inline void SetFactor(G4double val = 1.0)
  {
    if(writer != nullptr)
      writer->SetFactor(val);
  }
  inline G4double GetFactor() const
  {
    if(writer != nullptr)
    {
      return writer->GetFactor();
    }
    
    return -1.0;
  }

 protected:
  G4ScoringManager();

 private:
  // Disable copy constructor and assignement operator
  G4ScoringManager(const G4ScoringManager&);
  G4ScoringManager& operator=(const G4ScoringManager&);

 private:
  static G4ThreadLocal G4ScoringManager* fSManager;
  static G4ThreadLocal G4int replicaLevel;
  G4int verboseLevel;
  G4ScoringMessenger* fMessenger;
  G4ScoreQuantityMessenger* fQuantityMessenger;

  MeshVec fMeshVec;
  G4VScoringMesh* fCurrentMesh;

  G4VScoreWriter* writer;
  G4VScoreColorMap* fDefaultLinearColorMap;
  ColorMapDict* fColorMapDict;

  MeshMap fMeshMap;
};

#endif
