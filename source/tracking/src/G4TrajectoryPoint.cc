//
// ********************************************************************
// * DISCLAIMER                                                       *
// *                                                                  *
// * The following disclaimer summarizes all the specific disclaimers *
// * of contributors to this software. The specific disclaimers,which *
// * govern, are listed with their locations in:                      *
// *   http://cern.ch/geant4/license                                  *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.                                                             *
// *                                                                  *
// * This  code  implementation is the  intellectual property  of the *
// * GEANT4 collaboration.                                            *
// * By copying,  distributing  or modifying the Program (or any work *
// * based  on  the Program)  you indicate  your  acceptance of  this *
// * statement, and all its terms.                                    *
// ********************************************************************
//
//
// $Id: G4TrajectoryPoint.cc,v 1.15 2005-03-22 14:40:51 allison Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
//
// ---------------------------------------------------------------
//
// G4TrajectoryPoint.cc
//
// ---------------------------------------------------------------

#include "G4TrajectoryPoint.hh"

#include "G4AttDefStore.hh"
#include "G4AttDef.hh"
#include "G4AttValue.hh"
#include "G4UnitsTable.hh"

G4Allocator<G4TrajectoryPoint> aTrajectoryPointAllocator;

G4TrajectoryPoint::G4TrajectoryPoint()
{
  fPosition = G4ThreeVector(0.,0.,0.);
}

G4TrajectoryPoint::G4TrajectoryPoint(G4ThreeVector pos)
{
  fPosition = pos;
}

G4TrajectoryPoint::G4TrajectoryPoint(const G4TrajectoryPoint &right)
 : G4VTrajectoryPoint(),fPosition(right.fPosition)
{
}

G4TrajectoryPoint::~G4TrajectoryPoint()
{
}

const std::map<G4String,G4AttDef>* G4TrajectoryPoint::GetAttDefs() const
{
  G4bool isNew;
  std::map<G4String,G4AttDef>* store
    = G4AttDefStore::GetInstance("G4TrajectoryPoint",isNew);
  if (isNew) {
    G4String Pos("Pos");
    (*store)[Pos] = G4AttDef(Pos, "Position", "Physics","","G4ThreeVector");
  }
  return store;
}

std::vector<G4AttValue>* G4TrajectoryPoint::CreateAttValues() const
{
  std::vector<G4AttValue>* values = new std::vector<G4AttValue>;

  values->push_back(G4AttValue("Pos",G4BestUnit(fPosition,"Length"),""));

  return values;
}
