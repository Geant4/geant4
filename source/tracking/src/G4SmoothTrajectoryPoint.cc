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
// $Id: G4SmoothTrajectoryPoint.cc,v 1.8 2002-12-06 12:16:44 johna Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
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
#include "g4std/strstream"

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
						 G4std::vector<G4ThreeVector>* auxiliaryPoints)
: fPosition(pos),
  fAuxiliaryPointVector(auxiliaryPoints)
{}

G4SmoothTrajectoryPoint::G4SmoothTrajectoryPoint(const G4SmoothTrajectoryPoint &right)
: fPosition(right.fPosition),fAuxiliaryPointVector(right.fAuxiliaryPointVector)
{
}

G4SmoothTrajectoryPoint::~G4SmoothTrajectoryPoint()
{
  if(fAuxiliaryPointVector) {
    delete fAuxiliaryPointVector;
  }
}


const G4std::map<G4String,G4AttDef>*
G4SmoothTrajectoryPoint::GetAttDefs() const
{
  G4bool isNew;
  G4std::map<G4String,G4AttDef>* store
    = G4AttDefStore::GetInstance("G4SmoothTrajectoryPoint",isNew);
  if (isNew) {
    G4String Pos("Pos");
    (*store)[Pos] = G4AttDef(Pos, "Step Position",
			     "Physics","","G4ThreeVector");
    G4String Aux("Aux");
    (*store)[Aux] = G4AttDef(Aux, "Auxiliary Point Position",
			     "Physics","","G4ThreeVector");
  }
  return store;
}

G4std::vector<G4AttValue>* G4SmoothTrajectoryPoint::CreateAttValues() const
{
  char c[100];
  G4std::ostrstream s(c,100);

  G4std::vector<G4AttValue>* values = new G4std::vector<G4AttValue>;

  if (fAuxiliaryPointVector) {
    G4std::vector<G4ThreeVector>::iterator iAux;
    for (iAux = fAuxiliaryPointVector->begin();
	 iAux != fAuxiliaryPointVector->end(); ++iAux) {
      s.seekp(G4std::ios::beg);
      s << G4BestUnit(*iAux,"Length") << G4std::ends;
      values->push_back(G4AttValue("Aux",c,""));
    }
  }

  s.seekp(G4std::ios::beg);
  s << G4BestUnit(fPosition,"Length") << G4std::ends;
  values->push_back(G4AttValue("Pos",c,""));

  return values;
}
