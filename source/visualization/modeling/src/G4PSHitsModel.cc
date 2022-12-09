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
// Created:  Mar. 31, 2009  Akinori Kimura
// Model which draws the primtive scorer hits.
//

#include "G4PSHitsModel.hh"

#include "G4ModelingParameters.hh"
#include "G4VGraphicsScene.hh"
#include "G4Event.hh"
#include "G4ScoringManager.hh"
#include "G4VScoringMesh.hh"

G4PSHitsModel::~G4PSHitsModel () {}

G4PSHitsModel::G4PSHitsModel (const G4String& requestedMapName):
  fRequestedMapName(requestedMapName)
{
  fType = "G4PSHitsModel";
  fGlobalTag = "G4PSHitsModel for G4THitsMap<G4StatDouble> hits.";
  fGlobalDescription = fGlobalTag;
}

void G4PSHitsModel::DescribeYourselfTo (G4VGraphicsScene& sceneHandler)
{
  using MeshScoreMap = G4VScoringMesh::MeshScoreMap;
  using RunScore = G4VScoringMesh::RunScore;
  G4ScoringManager* scoringManager= G4ScoringManager::GetScoringManagerIfExist();
  if (scoringManager) {
    G4int nMeshes = (G4int)scoringManager->GetNumberOfMesh();
    for (G4int iMesh = 0; iMesh < nMeshes; ++iMesh) {
      G4VScoringMesh* mesh = scoringManager->GetMesh(iMesh);
      if (mesh && mesh->IsActive()) {
        MeshScoreMap scoreMap = mesh->GetScoreMap();
        for(MeshScoreMap::const_iterator i = scoreMap.cbegin();
            i != scoreMap.cend(); ++i) {
          const G4String& name = i->first;
          if (fRequestedMapName == "all" || name == fRequestedMapName) {
            RunScore* fpCurrentHits = i->second;
            if (fpCurrentHits) sceneHandler.AddCompound(*fpCurrentHits);
          }
        }
      }
    }
  }
}
