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
/// \file runAndEvent/RE04/src/RE04TrajectoryPoint.cc
/// \brief Implementation of the RE04TrajectoryPoint class
//
//
#include "RE04TrajectoryPoint.hh"

#include "G4AttDefStore.hh"
#include "G4AttDef.hh"
#include "G4AttValue.hh"
#include "G4UnitsTable.hh"

#include "G4Material.hh"

//#define G4ATTDEBUG
#ifdef G4ATTDEBUG
#include "G4AttCheck.hh"
#endif

G4ThreadLocal G4Allocator<RE04TrajectoryPoint> * faTrajPointAllocator = 0;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
RE04TrajectoryPoint::RE04TrajectoryPoint()
 : G4VTrajectoryPoint(),
   fPosition(G4ThreeVector(0.,0.,0.)), 
   fpMaterial(0)
{
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
RE04TrajectoryPoint::RE04TrajectoryPoint(G4ThreeVector pos,
                                         const G4Material* mat)
 : G4VTrajectoryPoint(),
   fPosition(pos), 
   fpMaterial(mat)
{
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
RE04TrajectoryPoint::RE04TrajectoryPoint(const RE04TrajectoryPoint &right)
 : G4VTrajectoryPoint(),
   fPosition(right.fPosition),
   fpMaterial(right.fpMaterial)
{
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
RE04TrajectoryPoint::~RE04TrajectoryPoint()
{
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
const std::map<G4String,G4AttDef>* RE04TrajectoryPoint::GetAttDefs() const
{
  G4bool isNew;
  std::map<G4String,G4AttDef>* store
    = G4AttDefStore::GetInstance("RE04TrajectoryPoint",isNew);
  if (isNew) {
    G4String pos("Pos");
    (*store)[pos] =
      G4AttDef(pos, "Position", "Physics","G4BestUnit","G4ThreeVector");
    G4String mat("Mat");
    (*store)[mat] =
      G4AttDef(mat, "Mat", "Physics","Name","G4String");
  }
  return store;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
std::vector<G4AttValue>* RE04TrajectoryPoint::CreateAttValues() const
{
  std::vector<G4AttValue>* values = new std::vector<G4AttValue>;

  values->push_back(G4AttValue("Pos",G4BestUnit(fPosition,"Length"),""));
  values->push_back(G4AttValue("Mat",fpMaterial->GetName(),""));

#ifdef G4ATTDEBUG
  G4cout << G4AttCheck(values,GetAttDefs());
#endif

  return values;
}
