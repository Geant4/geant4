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
/// \file ClusteringAlgo.hh
/// \brief Definition of the ClustreringAlgo class

#ifndef ClusteringAlgo_H
#define ClusteringAlgo_H 1

#include "ClusterSBPoints.hh"
#include "SBPoint.hh"

#include <map>

class ClusteringAlgoMessenger;

class ClusteringAlgo
{
public:

  ClusteringAlgo(G4double pEps, G4int pMinPts, G4double pSPointsProb,
      G4double pEMinDamage, G4double pEMaxDamage);
  ~ClusteringAlgo();

  // Get Set methods
  G4double GetEps()
  {
    return fEps;
  };
  void SetEps(G4double val)
  {
    fEps=val;
  };
  G4int GetMinPts()
  {
    return fMinPts;
  };
  void SetMinPts(G4int val)
  {
    fMinPts=val;
  };
  G4double GetSPointsProb()
  {
    return fSPointsProb;
  };
  void SetSPointsProb(G4double val)
  {
    fSPointsProb=val;
  };
  G4double GetEMinDamage()
  {
    return fEMinDamage;
  };
  void SetEMinDamage(G4double val)
  {
    fEMinDamage=val;
  };
  G4double GetEMaxDamage()
  {
    return fEMaxDamage;
  };
  void SetEMaxDamage(G4double val)
  {
    fEMaxDamage=val;
  };

  // Register a damage (position, edep)
  void RegisterDamage(G4ThreeVector, G4double);

  // Clustering Algorithm
  std::map<G4int,G4int> RunClustering();

  // Clean all data structures
  void  Purge();

  // Return the number of simple break
  G4int GetSSB() const;
  // Return the number of complex simple break
  G4int GetComplexSSB() const;
  // Return the number of double strand break
  G4int GetDSB() const;
  // Return a map representing cluster size distribution
  // first G4int : cluster size (1 = SSB)
  // second G4int : counts
  std::map<G4int,G4int> GetClusterSizeDistribution();

private:

  // Functions to check if SB candidate
  G4bool IsInSensitiveArea();
  G4bool IsEdepSufficient(G4double);
  // Check if a SB point can be merged to a cluster, and do it
  bool FindCluster(SBPoint* pPt);
  // Check if two points can be merged
  bool AreOnTheSameCluster(G4ThreeVector, G4ThreeVector,G4double);
  // Merge clusters
  void MergeClusters();
  // Add SSB to clusters
  void IncludeUnassociatedPoints();

  // Parameters to run clustering algorithm
  G4double fEps;         // distance to merge SBPoints
  G4int fMinPts;         // number of SBPoints to create a cluster
  G4double fSPointsProb; // probability for a point to be in the sensitive area
  G4double fEMinDamage;  // min energy to create a damage
  G4double fEMaxDamage;  // energy to have a probability to create a damage = 1

  // Data structure containing all SB points
  std::vector<SBPoint*> fpSetOfPoints;
  // Datya structure containing all clusters
  std::vector<ClusterSBPoints*> fpClusters;
  // ID of the next SB point
  unsigned int fNextSBPointID;

  ClusteringAlgoMessenger* fpClustAlgoMessenger;

};

#endif

