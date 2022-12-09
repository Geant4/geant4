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
// 
// John Allison  November 1996

#include "G4Polymarker.hh"

G4Polymarker::G4Polymarker ():
fMarkerType (G4Polymarker::dots)
{}

G4Polymarker::G4Polymarker (const G4VMarker& marker):
  G4VMarker (marker),
  fMarkerType (G4Polymarker::dots)
{}

G4Polymarker::~G4Polymarker () = default;

std::ostream& operator << (std::ostream& os, const G4Polymarker& marker) {
  os << "G4Polymarker: type: ";
  switch (marker.fMarkerType) {
  case G4Polymarker::dots:
    os << "dots"; break;
  case G4Polymarker::circles:
    os << "circles"; break;
  case G4Polymarker::squares:
    os << "squares"; break;
  default:
    os << "unrecognised"; break;
  }
  os << "\n  " << (G4VMarker) marker;
  os << "\n  " << (G4Point3DList) marker;
  return os;
}
