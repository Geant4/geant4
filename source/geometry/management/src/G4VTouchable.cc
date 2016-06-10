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
// $Id: G4VTouchable.cc 66356 2012-12-18 09:02:32Z gcosmo $
//
// 
// class G4VTouchable implementation
//
// --------------------------------------------------------------------

#include "G4VTouchable.hh"

G4VTouchable::G4VTouchable()
{
}

G4VTouchable::~G4VTouchable()
{
}

G4VPhysicalVolume* G4VTouchable::GetVolume(G4int) const
{
  G4Exception("G4VTouchable::GetVolume()", "GeomMgt0001",
              FatalException, "Undefined call to base class.");
  return 0;
}

G4VSolid* G4VTouchable::GetSolid(G4int) const
{
  G4Exception("G4VTouchable::GetSolid()", "GeomMgt0001",
              FatalException, "Undefined call to base class.");
  return 0;
}

G4int G4VTouchable::GetReplicaNumber(G4int) const
{
  G4Exception("G4VTouchable::GetReplicaNumber()", "GeomMgt0001",
              FatalException, "Undefined call to base class.");
  return 0;
}

G4int G4VTouchable::MoveUpHistory(G4int)
{
  G4Exception("G4VTouchable::MoveUpHistory()", "GeomMgt0001",
              FatalException, "Undefined call to base class.");
  return 0;
}

void G4VTouchable::UpdateYourself(G4VPhysicalVolume*,
			          const G4NavigationHistory* ) 
{
  G4Exception("G4VTouchable::UpdateYourself()", "GeomMgt0001",
              FatalException, "Undefined call to base class.");
}

G4int G4VTouchable::GetHistoryDepth() const
{
  G4Exception("G4VTouchable::GetHistoryDepth()", "GeomMgt0001",
              FatalException, "Undefined call to base class.");
  return  0;
}

const G4NavigationHistory* G4VTouchable::GetHistory() const
{
  G4Exception("G4VTouchable::GetHistory()", "GeomMgt0001",
              FatalException, "Undefined call to base class.");
  return 0;
}
