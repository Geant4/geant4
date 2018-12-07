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
/// \file ClusterSBPoints.cc
/// \brief Implementation of the ClusterSBPoints class

#include "ClusterSBPoints.hh"
#include "G4SystemOfUnits.hh"
#include <iostream>

using namespace std;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

ClusterSBPoints::ClusterSBPoints(std::set<SBPoint*> pSBPoints):
    fpRegisteredSBPoints()
{
  UpdateDoubleStrand();
  std::set<SBPoint*>::iterator itPt;
  for(itPt = pSBPoints.begin(); itPt != pSBPoints.end(); ++itPt)
  {
    AddSBPoint(*itPt);
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

ClusterSBPoints::~ClusterSBPoints()
{
  Clear();
}

void ClusterSBPoints::Clear()
{
  std::set<SBPoint*>::iterator itPt;
  for(itPt = fpRegisteredSBPoints.begin();
      itPt != fpRegisteredSBPoints.end();
      ++itPt)
  {
    (*itPt)->CleanCluster();
  }
  fpRegisteredSBPoints.clear();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void ClusterSBPoints::AddSBPoint(SBPoint* pSBPoint )
{
  assert(pSBPoint);
  fpRegisteredSBPoints.insert(pSBPoint);
  pSBPoint->SetCluster(this);

  UpdateDoubleStrand();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4ThreeVector ClusterSBPoints::GetBarycenter() const
{
  G4double x = 0;
  G4double y = 0;
  G4double z = 0;

  std::set<SBPoint* >::iterator itSDSPt;
  for(itSDSPt = fpRegisteredSBPoints.begin();
      itSDSPt != fpRegisteredSBPoints.end();
      ++itSDSPt)
  {
    x+= (*itSDSPt)->GetPosition().x();
    y+= (*itSDSPt)->GetPosition().y();
    z+= (*itSDSPt)->GetPosition().z();
  }

  return G4ThreeVector(
      x/fpRegisteredSBPoints.size(),
      y/fpRegisteredSBPoints.size(),
      z/fpRegisteredSBPoints.size());
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double ClusterSBPoints::GetEdep() const
{
  G4double res = 0;
  std::set<SBPoint* >::iterator itSDSPt;
  for(itSDSPt = fpRegisteredSBPoints.begin();
      itSDSPt != fpRegisteredSBPoints.end();
      ++itSDSPt)
  {
    res += (*itSDSPt)->GetEdep();
  }
  return res;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void ClusterSBPoints::UpdateDoubleStrand()
{
  fIsDoubleSB = false;
  bool firstStrandTouch = false;
  bool secondStrandTouch = false;

  std::set<SBPoint* >::iterator itSDSPt;
  for(itSDSPt = fpRegisteredSBPoints.begin();
      itSDSPt != fpRegisteredSBPoints.end();
      ++itSDSPt)
  {
    // if the SDSPoint is localized on the first strand
    if( ((*itSDSPt)->GetTouchedStrand() == 0 ) && !firstStrandTouch )
    {
      firstStrandTouch = true;
      if(secondStrandTouch)
      {
        fIsDoubleSB = true;
        return;
      }
    }
    // if the SDSPoint is localized on the second strand
    if( ((*itSDSPt)->GetTouchedStrand() == 1 ) && !secondStrandTouch )
    {
      secondStrandTouch = true;
      if(firstStrandTouch)
      {
        fIsDoubleSB = true;
        return;
      }
    }
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

bool AreOnTheSameCluster(const SBPoint* pPt1, const SBPoint* pPt2,
    G4double pMinDist)
{
  assert(pPt1);
  assert(pPt2);

  G4double x1=pPt1->GetPosition().x()/nm;
  G4double y1=pPt1->GetPosition().y()/nm;
  G4double z1=pPt1->GetPosition().z()/nm;

  G4double x2=pPt2->GetPosition().x()/nm;
  G4double y2=pPt2->GetPosition().y()/nm;
  G4double z2=pPt2->GetPosition().z()/nm;

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

void ClusterSBPoints::FindAllPointsPossible(std::vector<SBPoint*>* pPtsToCheck,
    G4int pMinPts, G4double pMinDist)
{
  assert((unsigned int)pMinPts > this->GetSize());
  std::vector<SBPoint*>::iterator itPt = pPtsToCheck->begin();
  while(itPt != pPtsToCheck->end() )
  {

    // If 1- each SBpoint is part of only one cluster
    //    2- the point isn't already in the cluster
    //    3- the point is close enough of the barycenter
    if( (!(*itPt)->HasCluster())
        && (fpRegisteredSBPoints.find(*itPt) == fpRegisteredSBPoints.end())
        && HasInBarycenter(*itPt, pMinDist)) // first version used HasIn method
    {
      // the point is added
      this->AddSBPoint(*itPt);
      if(this->GetSize() >= (unsigned int)pMinPts)
      {
        return;
      }
      // restart from scratch
      itPt = pPtsToCheck->begin();
    }else
    {
      ++itPt;
    }
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

bool ClusterSBPoints::HasIn(const SBPoint* pPtToCheck, G4double pMinDist)
{
  // check if the given point is near one of the cluster's point
  std::set<SBPoint*>::iterator itClusPt;
  for(itClusPt = fpRegisteredSBPoints.begin();
      itClusPt != fpRegisteredSBPoints.end();
      ++itClusPt)
  {
    // if are two different pts
    if((*pPtToCheck != *(*itClusPt)))
    {
      // if close enought of an include point of the cluster
      if( AreOnTheSameCluster(pPtToCheck, *itClusPt, pMinDist) )
      {
        return true;
      }
    }
  }
  return false;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

bool ClusterSBPoints::HasInBarycenter(const SBPoint* pPtToCheck,
    G4double pMinDist)
{

  G4double x1=pPtToCheck->GetPosition().x()/nm;
  G4double y1=pPtToCheck->GetPosition().y()/nm;
  G4double z1=pPtToCheck->GetPosition().z()/nm;

  G4double x2=this->GetBarycenter().x()/nm;
  G4double y2=this->GetBarycenter().y()/nm;
  G4double z2=this->GetBarycenter().z()/nm;

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

/// this will insert all registredSBPoint
/// from the given cluster to this cluster.
void ClusterSBPoints::MergeWith(ClusterSBPoints* pCluster)
{
  std::set<SBPoint*> points = pCluster->GetRegistredSBPoints();
  pCluster->Clear();
  std::set<SBPoint*>::iterator itPt;
  for(itPt = points.begin(); itPt != points.end(); ++itPt)
  {
    this->AddSBPoint(*itPt);
  }
}

