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
// $Id: G4Polyline.cc,v 1.9 2005/07/05 14:04:02 allison Exp $
// GEANT4 tag $Name: geant4-08-00 $
//
// 
// John Allison  July 1995

#include "G4Polyline.hh"

G4Polyline::G4Polyline () {}

G4Polyline::~G4Polyline () {}

G4Polyline& G4Polyline::transform (const G4Transform3D& transformation) {
  for (iterator i = begin(); i != end(); ++i) i->transform(transformation);
  return *this;
}

std::ostream& operator << (std::ostream& os, const G4Polyline& line) {
  os << "G4Polyline: ";
  os << '\n' << (G4Visible) line;
  os << '\n' << (G4Point3DList) line;
  return os;
}
