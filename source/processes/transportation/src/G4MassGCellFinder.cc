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
// $Id: G4MassGCellFinder.cc,v 1.3 2003/11/26 14:51:49 gcosmo Exp $
// GEANT4 tag $Name: geant4-08-00 $
//
// ----------------------------------------------------------------------
// GEANT 4 class source file
//
// G4MassGCellFinder.cc
//
// ----------------------------------------------------------------------

#include "G4MassGCellFinder.hh"
#include "G4TouchableHandle.hh"
#include "G4Step.hh"

G4MassGCellFinder::G4MassGCellFinder()
{
}

G4MassGCellFinder::~G4MassGCellFinder()
{
}

G4GeometryCell
G4MassGCellFinder::GetPreGeometryCell(const G4Step &aStep) const
{
  const G4TouchableHandle& th =
    aStep.GetPreStepPoint()->GetTouchableHandle();
  G4GeometryCell g(*th->GetVolume(), th->GetReplicaNumber()); 
  return g;
}

G4GeometryCell
G4MassGCellFinder::GetPostGeometryCell(const G4Step &aStep) const
{
  const G4TouchableHandle& th =
    aStep.GetPostStepPoint()->GetTouchableHandle();
  G4GeometryCell g(*th->GetVolume(), th->GetReplicaNumber());    
  return g;
}
