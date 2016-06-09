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
//
// $Id: G4MassGCellFinder.cc,v 1.5 2006/11/14 09:11:18 gcosmo Exp $
// GEANT4 tag $Name: geant4-08-02 $
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
