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
// $Id: G4TrajectoryPoint.cc,v 1.8 2002-10-16 11:38:38 johna Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
//
// ---------------------------------------------------------------
//
// G4TrajectoryPoint.cc
//
// ---------------------------------------------------------------

#include "G4TrajectoryPoint.hh"

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
 : fPosition(right.fPosition)
{
}

G4TrajectoryPoint::~G4TrajectoryPoint()
{
}

const G4std::vector<G4AttDef>* G4TrajectoryPoint::GetAttDefs() const
{ return 0; }  // Empty for now.

G4std::vector<G4AttValue>* G4TrajectoryPoint::GetAttValues() const
{ return 0; }  // Empty for now.
