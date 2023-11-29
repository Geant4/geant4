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
// This example is provided by the Geant4-DNA collaboration
// Any report or published results obtained using the Geant4-DNA software
// shall cite the following Geant4-DNA collaboration publication:
// Med. Phys. 37 (2010) 4692-4708
// The Geant4-DNA web site is available at http://geant4-dna.org
//
// Authors: Henri Payno and Yann Perrot
//
//
/// \file ClusteringAlgo.cc
/// \brief Implementation of the ClustreringAlgo class

#include "ClusteringAlgo.hh"
#include "ClusteringAlgoMessenger.hh"

#include "G4SystemOfUnits.hh"
#include "Randomize.hh"

#include <map>

using namespace std;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

ClusteringAlgo::ClusteringAlgo(G4double pEps,G4int pMinPts,
    G4double pSPointsProb,G4double pEMinDamage, G4double pEMaxDamage)
:fEps(pEps),fMinPts(pMinPts),
 fSPointsProb(pSPointsProb),fEMinDamage(pEMinDamage),fEMaxDamage(pEMaxDamage)
{
  fNextSBPointID = 0;
  fpClustAlgoMessenger = new ClusteringAlgoMessenger(this);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

ClusteringAlgo::~ClusteringAlgo()
{
  delete fpClustAlgoMessenger;
  Purge();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

// Random sampling in space
G4bool ClusteringAlgo::IsInSensitiveArea()
{
  return fSPointsProb > G4UniformRand();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

// Random sampling in energy
G4bool ClusteringAlgo::IsEdepSufficient(G4double pEdep)
{
  if(pEdep<fEMinDamage)
  {
    return false;
  }

  else if(pEdep>fEMaxDamage)
  {
    return true;
  }
  else
  {
    G4double proba = (pEdep/eV - fEMinDamage/eV)/
        (fEMaxDamage/eV-fEMinDamage/eV);
    return (proba>G4UniformRand());
  }

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

// Add an event interaction to the unregistered damage if
// good conditions (pos and energy) are met
//

void ClusteringAlgo::RegisterDamage(G4ThreeVector pPos,G4double pEdep)
{
  if(IsEdepSufficient(pEdep))
  {
    if(IsInSensitiveArea())
    {
      fpSetOfPoints.push_back( new SBPoint(fNextSBPointID++, pPos, pEdep));
    }
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

map<G4int,G4int> ClusteringAlgo::RunClustering()
{

  // quick sort style
  // create cluster
  std::vector<SBPoint*>::iterator itVisitorPt, itObservedPt;
  for(itVisitorPt = fpSetOfPoints.begin();
      itVisitorPt != fpSetOfPoints.end();
      ++itVisitorPt  )
  {
    itObservedPt = itVisitorPt;
    itObservedPt ++;
    while(itObservedPt != fpSetOfPoints.end() )
    {
      // if at least one of the two points has not a cluster
      if(!((*itObservedPt)->HasCluster() && (*itVisitorPt)->HasCluster()))
      {
        if(AreOnTheSameCluster( (*itObservedPt)->GetPosition(),
                                (*itVisitorPt)->GetPosition(),fEps))
        {
          // if none has a cluster. Create a new one
          if(!(*itObservedPt)->HasCluster() && !(*itVisitorPt)->HasCluster())
          {
            // create the new cluster
            set<SBPoint*> clusterPoints;
            clusterPoints.insert((*itObservedPt));
            clusterPoints.insert((*itVisitorPt));
            ClusterSBPoints* lCluster = new ClusterSBPoints(clusterPoints);
            assert(lCluster);
            fpClusters.push_back(lCluster);
            assert(lCluster);
            // inform SB point that they are part of a cluster now
            assert(lCluster);
            (*itObservedPt)->SetCluster(lCluster);
            assert(lCluster);
            (*itVisitorPt)->SetCluster(lCluster);
          }else
          {
            // add the point to the existing cluster
            if((*itObservedPt)->HasCluster())
            {
              (*itObservedPt)->GetCluster()->AddSBPoint((*itVisitorPt));
              (*itVisitorPt)->SetCluster((*itObservedPt)->GetCluster());
            }

            if((*itVisitorPt)->HasCluster())
            {
              (*itVisitorPt)->GetCluster()->AddSBPoint((*itObservedPt));
              (*itObservedPt)->SetCluster((*itVisitorPt)->GetCluster());
            }
          }
        }
      }
      ++itObservedPt;
    }
  }

  // associate isolated points and merge clusters
  IncludeUnassociatedPoints();
  MergeClusters();

  // return cluster size distribution
  return GetClusterSizeDistribution();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

// try to merge cluster between them, based on the distance between barycenters
void ClusteringAlgo::MergeClusters()
{
  std::vector<ClusterSBPoints*>::iterator itCluster1, itCluster2;
  for(itCluster1 = fpClusters.begin();
      itCluster1 != fpClusters.end();
      ++itCluster1)
  {
    G4ThreeVector baryCenterClust1 = (*itCluster1)->GetBarycenter();
    itCluster2 = itCluster1;
    itCluster2++;
    while(itCluster2 != fpClusters.end())
    {
      G4ThreeVector baryCenterClust2 = (*itCluster2)->GetBarycenter();
      // if we can merge both cluster
      if(AreOnTheSameCluster(baryCenterClust1, baryCenterClust2,fEps))
      {
        (*itCluster1)->MergeWith(*itCluster2);
        delete *itCluster2;
        fpClusters.erase(itCluster2);
        return MergeClusters();
      }else
      {
        itCluster2++;
      }
    }
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void ClusteringAlgo::IncludeUnassociatedPoints()
{
  std::vector<SBPoint*>::iterator itVisitorPt;
  // Associate all point not in a cluster if possible ( to the first found cluster)
  for(itVisitorPt = fpSetOfPoints.begin();
      itVisitorPt != fpSetOfPoints.end();
      ++itVisitorPt)
  {
    if(!(*itVisitorPt)->HasCluster())
    {
      FindCluster(*itVisitorPt);
    }
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

bool ClusteringAlgo::FindCluster(SBPoint* pPt)
{
  assert(!pPt->HasCluster());
  std::vector<ClusterSBPoints*>::iterator itCluster;
  for(itCluster = fpClusters.begin();
      itCluster != fpClusters.end();
      ++itCluster)
  {
    //if((*itCluster)->hasIn(pPt, fEps))
    if((*itCluster)->HasInBarycenter(pPt, fEps))
    {
      (*itCluster)->AddSBPoint(pPt);
      return true;
    }
  }  
  return false;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

bool ClusteringAlgo::AreOnTheSameCluster(G4ThreeVector pPt1,
    G4ThreeVector pPt2, G4double pMinDist)
{
  G4double x1=pPt1.x()/nm;
  G4double y1=pPt1.y()/nm;
  G4double z1=pPt1.z()/nm;

  G4double x2=pPt2.x()/nm;
  G4double y2=pPt2.y()/nm;
  G4double z2=pPt2.z()/nm;

  // if the two points are closed enough
  if(((x1-x2)*(x1-x2)+(y1-y2)*(y1-y2)+(z1-z2)*(z1-z2))<=
      (pMinDist/nm*pMinDist/nm))
  {
    return true;
  }else
  {
    return false;
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4int ClusteringAlgo::GetSSB() const
{
  G4int nbSSB = 0;
  std::vector<SBPoint*>::const_iterator itSDSPt;
  for(itSDSPt = fpSetOfPoints.begin();
      itSDSPt != fpSetOfPoints.end();
      ++itSDSPt)
  {
    if(!(*itSDSPt)->HasCluster())
    {
      nbSSB++;
    }
  }
  return nbSSB;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4int ClusteringAlgo::GetComplexSSB() const
{
  G4int nbSSB = 0;
  std::vector<ClusterSBPoints*>::const_iterator itCluster;
  for(itCluster = fpClusters.begin();
      itCluster != fpClusters.end();
      ++itCluster)
  {
    if((*itCluster)->IsSSB())
    {
      nbSSB ++;
    }
  }
  return nbSSB;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4int ClusteringAlgo::GetDSB() const
{
  G4int nbDSB = 0;
  std::vector<ClusterSBPoints*>::const_iterator itCluster;
  for(itCluster = fpClusters.begin();
      itCluster != fpClusters.end();
      ++itCluster)
  {
    if((*itCluster)->IsDSB())
    {
      nbDSB ++;
    }
  }
  return nbDSB;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

map<G4int,G4int> ClusteringAlgo::GetClusterSizeDistribution()
{
  std::map<G4int,G4int>  sizeDistribution;
  sizeDistribution[1]=GetSSB();
  std::vector<ClusterSBPoints*>::const_iterator itCluster;
  for(itCluster = fpClusters.begin();
      itCluster != fpClusters.end();
      itCluster++)
  {
    sizeDistribution[(*itCluster)->GetSize()]++;
  }
  return sizeDistribution;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void ClusteringAlgo::Purge()
{
  fNextSBPointID = 0;
  std::vector<ClusterSBPoints*>::iterator itCluster;
  for(itCluster = fpClusters.begin();
      itCluster != fpClusters.end();
      ++itCluster)
  {
    delete *itCluster;
    *itCluster = NULL;
  }  
  fpClusters.clear();
  std::vector<SBPoint*>::iterator itPt;
  for(itPt = fpSetOfPoints.begin();
      itPt != fpSetOfPoints.end();
      ++itPt)
  {
    delete *itPt;
    *itPt = NULL;
  }
  fpSetOfPoints.clear();
}

