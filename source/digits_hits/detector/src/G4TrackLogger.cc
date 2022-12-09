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
//
// ----------------------------------------------------------------------
// GEANT 4 class source file
//
// G4TrackLogger.cc
//
// ----------------------------------------------------------------------

#include "G4TrackLogger.hh"

G4TrackLogger::G4TrackLogger()
  : fPreviousEventID(-1)
{}

G4TrackLogger::~G4TrackLogger() {}

void G4TrackLogger::SetEventID(G4int id)
{
  if(id != fPreviousEventID)
  {
    fTrackIDsSet.clear();
    fPreviousEventID = id;
  }
}

G4bool G4TrackLogger::FirstEnterance(G4int trid)
{
  G4bool first = true;
  auto n       = fTrackIDsSet.count(trid);
  if(n == 1)
  {
    first = false;
  }
  else if(n == 0)
  {
    fTrackIDsSet.insert(trid);
  }
  else if(n > 1)
  {
    G4cout << "Error G4TrackLogger::FirstEnterance: "
           << "more than one elm in set!" << G4endl;
  }
  return first;
}
