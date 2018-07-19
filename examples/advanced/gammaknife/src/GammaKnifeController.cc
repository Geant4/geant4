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


#include "GammaKnifeController.hh"
#include "G4UImanager.hh"
#include "G4RunManager.hh"
#include <fstream>

#include "G4SystemOfUnits.hh"

GammaKnifeController::GammaKnifeController( GammaKnifeDetectorConstruction *det ) : scoreMaps( 0 )
{
    detector = det;

    messenger = new GammaKnifeMessenger( this );
}

GammaKnifeController::~GammaKnifeController()
{
    delete messenger;
}

void GammaKnifeController::BeamOn( G4int n_event )
{
    PrepareHitsAccumulation();
    for (G4int i = 0; i < GAMMAKNIFE_SOURCES; i++)
    {
        RotateForward(i);
	G4RunManager::GetRunManager()->BeamOn(n_event);

        if (i != (GAMMAKNIFE_SOURCES - 1))        
	  StoreHits();
        
        RotateBack(i);
    }
    AccumulateAllHits();
}

void GammaKnifeController::RotateForward( G4int position )
{
    // Rotate all scoring meshes to the right position
    G4ScoringManager* scm = G4ScoringManager::GetScoringManagerIfExist();
    if (scm)
    {
        for (size_t i = 0; i < scm->GetNumberOfMesh(); i++)
        {
            G4VScoringMesh * mesh = scm->GetMesh(i);
            mesh->RotateX( thetaAngles[position] );
            mesh->RotateZ( phiAngles[position] );
        }
    }
}

void GammaKnifeController::RotateBack( G4int position )
{
    G4ScoringManager* scm = G4ScoringManager::GetScoringManagerIfExist();
    if (scm)
    {
        for (size_t i = 0; i < scm->GetNumberOfMesh(); i++)
        {
            G4VScoringMesh * mesh = scm->GetMesh(i);
            mesh->RotateZ( - phiAngles[position] );
            mesh->RotateX( - thetaAngles[position] );
        }
    }
}

void GammaKnifeController::PrepareHitsAccumulation()
{
    G4ScoringManager* scm = G4ScoringManager::GetScoringManagerIfExist();
    if (scm)
    {
        size_t size = scm->GetNumberOfMesh();
        scoreMaps = new MeshScoreMap[size];
        for (size_t i = 0; i < size; i++)
        {
            G4VScoringMesh * mesh = scm->GetMesh(i);

            MeshScoreMap scoreMap = mesh->GetScoreMap();
            MeshScoreMap& storedScoreMap = scoreMaps[i];

            MeshScoreMap::iterator it = scoreMap.begin();
            for( ; it != scoreMap.end(); it++)
            {
                std::string hitMapName = it->first;
                G4THitsMap<G4StatDouble>* hitMapToStore
                   = new G4THitsMap<G4StatDouble>("GammaKnifeController", hitMapName);
                storedScoreMap[ hitMapName ] = hitMapToStore;
            }
        }
    }
}

void GammaKnifeController::StoreHits()
{
    G4ScoringManager* scm = G4ScoringManager::GetScoringManagerIfExist();
    if (scm)
    {
        for (size_t i = 0; i < scm->GetNumberOfMesh(); i++)
        {
            G4VScoringMesh* mesh = scm->GetMesh(i);


            MeshScoreMap scoreMap = mesh->GetScoreMap();
            MeshScoreMap& storedScoreMap = scoreMaps[i];

            MeshScoreMap::iterator it = scoreMap.begin();
            for( ; it != scoreMap.end(); it++)
            {
               std::string hitMapName = it->first;
               //*storedScoreMap[hitMapName] += *(it->second);
               auto storedMap = storedScoreMap[hitMapName]->GetMap();
               auto mapItr = it->second->GetMap()->begin();
               for(;mapItr!=it->second->GetMap()->end();mapItr++)
               {
                 auto key = mapItr->first;
                 auto val = mapItr->second;
                 if(storedMap->find(key)==storedMap->end())
                 { (*storedMap)[key] = new G4StatDouble(); }
                 (*storedMap)[key]->add(val);
               }
            }
        }
    }
}

void GammaKnifeController::AccumulateAllHits()
{
    G4ScoringManager* scm = G4ScoringManager::GetScoringManagerIfExist();
    if (scm)
    {
        for (size_t i = 0; i < scm->GetNumberOfMesh(); i++)
        {
            G4VScoringMesh* mesh = scm->GetMesh(i);
            MeshScoreMap& storedScoreMap = scoreMaps[i];
            MeshScoreMap::iterator it = storedScoreMap.begin();
            for( ; it != storedScoreMap.end(); it++)
            {
                mesh->Accumulate( it->second );
            }
        }
    }
}

void GammaKnifeController::ReadFile( std::string fileName )
{
  //G4cout << "Enter ReadFile()...";
  const int SZ = 100;
  char buf[SZ];
  
  phiAngles.clear();    // If called for the second time
  thetaAngles.clear();  // we won't have 402 positions...
  
  std::ifstream ifs;
  ifs.open( fileName.c_str() );
  
  for (G4int i = 0; i < GAMMAKNIFE_SOURCES; i++)
    {
      G4double phi, theta;
      
      /* Skip the "Axx" at the beginning of the line */
      for (G4int c = 0; c < 4; c++) ifs.get();
      
      ifs >> phi >> theta;
      ifs.getline(buf, SZ); // Next line
      
      phiAngles.push_back( phi * degree );
      thetaAngles.push_back( theta * degree );
    }
  ifs.close();
  //G4cout << "... done " << G4endl;
}
