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
// $Id: G4ParallelGCellFinder.cc,v 1.1 2002-10-16 16:27:41 dressel Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// ----------------------------------------------------------------------
// GEANT 4 class source file
//
// G4ParallelGCellFinder.cc
//
// ----------------------------------------------------------------------

#include "G4ParallelGCellFinder.hh"
#include "G4VParallelStepper.hh"


G4ParallelGCellFinder::
G4ParallelGCellFinder(const G4VParallelStepper &astepper):
  fPStepper(astepper)
{}

G4ParallelGCellFinder::~G4ParallelGCellFinder()
{}

G4GeometryCell G4ParallelGCellFinder::
GetPreGeometryCell(const G4Step &aStep) const {
  G4GeometryCell g(fPStepper.GetPStep().GetPreGeometryCell());
  return g;
}

G4GeometryCell G4ParallelGCellFinder::
GetPostGeometryCell(const G4Step &aStep) const {
  G4GeometryCell g(fPStepper.GetPStep().GetPostGeometryCell());
  return g;
}
