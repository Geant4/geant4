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
// G4ClonedSmoothTrajectoryPoint class implementation
//
// Makoto Asai (JLab) - Oct.2024
// --------------------------------------------------------------------

#include "G4ClonedSmoothTrajectoryPoint.hh"
#include "G4SmoothTrajectoryPoint.hh"

#include "G4AttDef.hh"
#include "G4AttDefStore.hh"
#include "G4AttValue.hh"
#include "G4UnitsTable.hh"

// #define G4ATTDEBUG
#ifdef G4ATTDEBUG
#  include "G4AttCheck.hh"
#endif

G4Allocator<G4ClonedSmoothTrajectoryPoint>*& aClonedSmoothTrajectoryPointAllocator()
{
  static G4Allocator<G4ClonedSmoothTrajectoryPoint>* _instance = nullptr;
  return _instance;
}

G4ClonedSmoothTrajectoryPoint::G4ClonedSmoothTrajectoryPoint(
  G4ThreeVector pos, std::vector<G4ThreeVector>* auxiliaryPoints)
  : fPosition(pos), fAuxiliaryPointVector(auxiliaryPoints)
{}

G4ClonedSmoothTrajectoryPoint::G4ClonedSmoothTrajectoryPoint(const G4SmoothTrajectoryPoint& right)
  : fPosition(right.fPosition), fAuxiliaryPointVector(right.fAuxiliaryPointVector)
{}

G4ClonedSmoothTrajectoryPoint::~G4ClonedSmoothTrajectoryPoint() { delete fAuxiliaryPointVector; }

const std::map<G4String, G4AttDef>* G4ClonedSmoothTrajectoryPoint::GetAttDefs() const
{
  G4bool isNew;
  std::map<G4String, G4AttDef>* store =
    G4AttDefStore::GetInstance("G4ClonedSmoothTrajectoryPoint", isNew);
  if (isNew) {
    G4String Pos("Pos");
    (*store)[Pos] = G4AttDef(Pos, "Step Position", "Physics", "G4BestUnit", "G4ThreeVector");
    G4String Aux("Aux");
    (*store)[Aux] =
      G4AttDef(Aux, "Auxiliary Point Position", "Physics", "G4BestUnit", "G4ThreeVector");
  }
  return store;
}

std::vector<G4AttValue>* G4ClonedSmoothTrajectoryPoint::CreateAttValues() const
{
  auto values = new std::vector<G4AttValue>;

  if (fAuxiliaryPointVector != nullptr) {
    for (const auto& iAux : *fAuxiliaryPointVector) {
      values->push_back(G4AttValue("Aux", G4BestUnit(iAux, "Length"), ""));
    }
  }

  values->push_back(G4AttValue("Pos", G4BestUnit(fPosition, "Length"), ""));

#ifdef G4ATTDEBUG
  G4cout << G4AttCheck(values, GetAttDefs());
#endif

  return values;
}
