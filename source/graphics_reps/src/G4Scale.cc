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
// $Id: G4Scale.cc,v 1.3 2001-07-24 21:41:43 johna Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// John Allison  21st July 2001

#include "G4Scale.hh"

G4Scale::G4Scale (G4double length, const G4String& annotation,
		  G4Scale::Direction direction,
		  G4bool autoPlacing,
		  G4double xmid, G4double ymid, G4double zmid):
  fLength(length),
  fAnnotation(annotation),
  fDirection(direction),
  fAutoPlacing(autoPlacing),
  fXmid(xmid), fYmid(ymid), fZmid(zmid) {}

G4Scale::~G4Scale () {}

G4Visible & G4Scale::operator = (const G4Visible &from) {
  return G4Visible::operator = (from);
}

G4VVisPrim & G4Scale::operator = (const G4VVisPrim &from) {
  return G4VVisPrim::operator = (from);
}

G4VMarker & G4Scale::operator = (const G4VMarker &from) {
  return G4VMarker::operator = (from);
}

G4Scale & G4Scale::operator = (const G4Scale &from) {
  if (&from == this) return *this;
  G4VMarker::operator = (from);
  fLength = from.fLength;
  fAnnotation = from.fAnnotation;
  fDirection = from.fDirection;
  fAutoPlacing = from.fAutoPlacing;
  fXmid = from.fXmid;
  fYmid = from.fYmid;
  fZmid = from.fZmid;
  return *this;
}
