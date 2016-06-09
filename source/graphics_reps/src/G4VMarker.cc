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
// $Id: G4VMarker.cc,v 1.9 2005/07/05 14:04:02 allison Exp $
// GEANT4 tag $Name: geant4-08-00 $
//
// 

#include "G4VMarker.hh"

#include "G4VisAttributes.hh"

G4VMarker::G4VMarker ():
  fPosition   (G4Point3D ()),
  fWorldSize  (0.),
  fScreenSize (0.),
  fFillStyle  (noFill),
  fInfo       (G4String())
{}

G4VMarker::G4VMarker (const G4Point3D& pos):
  fPosition   (pos),
  fWorldSize  (0.),
  fScreenSize (0.),
  fFillStyle  (noFill),
  fInfo       (G4String())
{}

G4VMarker::~G4VMarker () {}

G4bool G4VMarker::operator != (const G4VMarker& m) const {
  if (
      (G4Visible::operator != (m))   ||
      (fWorldSize  != m.fWorldSize)   ||
      (fScreenSize != m.fScreenSize)  ||
      (fFillStyle  != m.fFillStyle)   ||
      !(fPosition  == m.fPosition)
      )
    return true;
  return false;
}

std::ostream& operator << (std::ostream& os, const G4VMarker& marker) {
  os << "G4VMarker: position: " << marker.fPosition
     << ", world size: " << marker.fWorldSize
     << ", screen size: " << marker.fScreenSize << '\n'
     << "           fill style: ";
  switch (marker.fFillStyle) {
  case G4VMarker::noFill:
    os << "no fill";
    break;
  case G4VMarker::hashed:
    os << "hashed";
    break;
  case G4VMarker::filled:
    os << "filled";
    break;
  default:
    os << "unrecognised"; break;
  }
  os << "\n           " << (G4Visible) marker;
  return os;
}

const G4String& G4VMarker::GetInfo() const { return fInfo ;}

void G4VMarker::SetInfo( const G4String& info ){ fInfo = info ;}

