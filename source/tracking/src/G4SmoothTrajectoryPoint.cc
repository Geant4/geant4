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
// $Id: G4SmoothTrajectoryPoint.cc,v 1.1 2002-09-04 02:09:38 asaim Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
//
// ---------------------------------------------------------------
//
// G4SmoothTrajectoryPoint.cc
//
// ---------------------------------------------------------------

#include "G4SmoothTrajectoryPoint.hh"

G4Allocator<G4SmoothTrajectoryPoint> aTrajectoryPointAllocator;

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

G4SmoothTrajectoryPoint::G4SmoothTrajectoryPoint(const G4SmoothTrajectoryPoint &right)
: fPosition(right.fPosition),fAuxiliaryPointVector(right.fAuxiliaryPointVector)
{
}

G4SmoothTrajectoryPoint::~G4SmoothTrajectoryPoint()
{
  if(!fAuxiliaryPointVector) {
    for(size_t i=0;i<fAuxiliaryPointVector->size();i++) {
     delete (*fAuxiliaryPointVector)[i];
    }
    fAuxiliaryPointVector->clear();
    delete fAuxiliaryPointVector;
  }
}

const G4AttValueList* G4SmoothTrajectoryPoint::GetAttValues() const
{ return 0; }

