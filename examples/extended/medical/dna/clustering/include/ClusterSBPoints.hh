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
/// \file ClusterSBPoints.hh
/// \brief Definition of the ClusterSBPoints class

#ifndef CLUSTER_SB_POINT_HH
#define CLUSTER_SB_POINT_HH

#include "SBPoint.hh"

#include <set>
#include <vector>

/// \brief define a cluster of SB Points
class ClusterSBPoints
{
public:
  ClusterSBPoints(std::set<SBPoint*> pPoints);
  ~ClusterSBPoints();

  bool IsDSB() const
  {
    return fIsDoubleSB;
  }

  bool IsSSB() const
  {
    return ! IsDSB();
  }

  unsigned int GetSize() const
  {
    return fpRegisteredSBPoints.size();
  }

  void AddSBPoint(SBPoint* pSBPoint );
  G4ThreeVector GetBarycenter() const;

  G4double GetEdep() const;

  void FindAllPointsPossible(std::vector<SBPoint*>* ptToCheck,
      G4int minPts, G4double minDist);

  /// \brief will set to all store SBPoint they are part of this cluster
  void NoticeInternalPts();

  bool HasIn(const SBPoint*, G4double);
  bool HasInBarycenter(const SBPoint*, G4double);

  void MergeWith(ClusterSBPoints*);

  std::set<SBPoint* > GetRegistredSBPoints() const
      {
    return fpRegisteredSBPoints;
      }

  void Clear();

private:
  void UpdateDoubleStrand();

  unsigned int fSize;
  bool fIsDoubleSB; // is a double strand break ?

  std::set<SBPoint* > fpRegisteredSBPoints;
};

#endif //CLUSTER_SB_POINT_HH
