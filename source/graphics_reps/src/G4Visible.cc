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
// $Id: G4Visible.cc,v 1.13 2005/07/05 14:04:02 allison Exp $
// GEANT4 tag $Name: geant4-08-00 $
//
// 
// John Allison  30th October 1996
// Base class for all things visible, i.e., which have Vis Attributes.

#include "G4Visible.hh"
#include "G4VisAttributes.hh"
#include "G4ios.hh"

G4Visible::G4Visible (): fpVisAttributes (0) {}

G4Visible::~G4Visible () {}

G4Visible::G4Visible (const G4VisAttributes* pVA):
  fpVisAttributes (pVA)
{}

void G4Visible::SetVisAttributes (const G4VisAttributes& VA) {
  fpVisAttributes = new G4VisAttributes(VA);
}

G4bool G4Visible::operator != (const G4Visible& right) const {
  // Simple test on non-equality of address...
  return fpVisAttributes != right.fpVisAttributes;
}

std::ostream& operator << (std::ostream& os, const G4Visible& v) {
  if (v.fpVisAttributes) return os << *(v.fpVisAttributes);
  else return os << "No Visualization Attributes";
}
