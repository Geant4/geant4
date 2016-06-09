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
// $Id: G4SmoothTrajectoryPoint.cc,v 1.14 2005/05/03 17:48:51 allison Exp $
// GEANT4 tag $Name: geant4-08-00 $
//
//
// ---------------------------------------------------------------
//
// G4SmoothTrajectoryPoint.cc
//
// ---------------------------------------------------------------

#include "G4SmoothTrajectoryPoint.hh"

#include "G4AttDefStore.hh"
#include "G4AttDef.hh"
#include "G4AttValue.hh"
#include "G4UnitsTable.hh"

//#define G4ATTDEBUG
#ifdef G4ATTDEBUG
#include "G4AttCheck.hh"
#endif

G4Allocator<G4SmoothTrajectoryPoint> aSmoothTrajectoryPointAllocator;

G4SmoothTrajectoryPoint::G4SmoothTrajectoryPoint()
: fAuxiliaryPointVector(0)
{
  fPosition = G4ThreeVector(0.,0.,0.);
}

G4SmoothTrajectoryPoint::G4SmoothTrajectoryPoint(G4ThreeVector pos)
: fAuxiliaryPointVector(0)
{
  fPosition = pos;
}

G4SmoothTrajectoryPoint::G4SmoothTrajectoryPoint(G4ThreeVector pos,
						 std::vector<G4ThreeVector>* auxiliaryPoints)
: fPosition(pos),
  fAuxiliaryPointVector(auxiliaryPoints)
{}

G4SmoothTrajectoryPoint::G4SmoothTrajectoryPoint(const G4SmoothTrajectoryPoint &right)
: G4VTrajectoryPoint(),
  fPosition(right.fPosition),fAuxiliaryPointVector(right.fAuxiliaryPointVector)
{
}

G4SmoothTrajectoryPoint::~G4SmoothTrajectoryPoint()
{
  if(fAuxiliaryPointVector) {
    delete fAuxiliaryPointVector;
  }
}


const std::map<G4String,G4AttDef>*
G4SmoothTrajectoryPoint::GetAttDefs() const
{
  G4bool isNew;
  std::map<G4String,G4AttDef>* store
    = G4AttDefStore::GetInstance("G4SmoothTrajectoryPoint",isNew);
  if (isNew) {
    G4String Pos("Pos");
    (*store)[Pos] = G4AttDef(Pos, "Step Position",
			     "Physics","G4BestUnit","G4ThreeVector");
    G4String Aux("Aux");
    (*store)[Aux] = G4AttDef(Aux, "Auxiliary Point Position",
			     "Physics","G4BestUnit","G4ThreeVector");
  }
  return store;
}

std::vector<G4AttValue>* G4SmoothTrajectoryPoint::CreateAttValues() const
{
  std::vector<G4AttValue>* values = new std::vector<G4AttValue>;

  if (fAuxiliaryPointVector) {
    std::vector<G4ThreeVector>::iterator iAux;
    for (iAux = fAuxiliaryPointVector->begin();
	 iAux != fAuxiliaryPointVector->end(); ++iAux) {
      values->push_back(G4AttValue("Aux",G4BestUnit(*iAux,"Length"),""));
    }
  }

  values->push_back(G4AttValue("Pos",G4BestUnit(fPosition,"Length"),""));

#ifdef G4ATTDEBUG
  G4cout << G4AttCheck(values,GetAttDefs());
#endif

  return values;
}
