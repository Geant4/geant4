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
// John Allison  30th October 1996
// Base class for all things visible, i.e., which have Vis Attributes.

#include "G4Visible.hh"
#include "G4VisAttributes.hh"
#include "G4ios.hh"

G4Visible::G4Visible ():
  fpVisAttributes (nullptr),
  fAllocatedVisAttributes (false)
{}

G4Visible::G4Visible (const G4Visible& visible){
  fAllocatedVisAttributes = visible.fAllocatedVisAttributes;
  if (fAllocatedVisAttributes)
    fpVisAttributes = new G4VisAttributes(*visible.fpVisAttributes);
  else fpVisAttributes = visible.fpVisAttributes;
}

G4Visible::G4Visible (G4Visible&& visible){
  fAllocatedVisAttributes = visible.fAllocatedVisAttributes;
  fpVisAttributes = visible.fpVisAttributes;
  visible.fpVisAttributes = nullptr;
  visible.fAllocatedVisAttributes = false;
}

G4Visible::G4Visible (const G4VisAttributes* pVA):
  fpVisAttributes (pVA),
  fAllocatedVisAttributes (false)
{}

G4Visible::~G4Visible () {
  if (fAllocatedVisAttributes) delete fpVisAttributes;
}

G4Visible& G4Visible::operator= (const G4Visible& rhs) {
  if (&rhs == this) return *this;
  fInfo = rhs.fInfo;
  fAllocatedVisAttributes = rhs.fAllocatedVisAttributes;
  if (fAllocatedVisAttributes) {
    delete fpVisAttributes;
    fpVisAttributes = new G4VisAttributes(*rhs.fpVisAttributes);
  }
  else fpVisAttributes = rhs.fpVisAttributes;
  return *this;
}

G4Visible& G4Visible::operator= (G4Visible&& rhs) {
  if (&rhs == this) return *this;
  fInfo = rhs.fInfo;
  if (fAllocatedVisAttributes) delete fpVisAttributes;
  fpVisAttributes = rhs.fpVisAttributes;
  fAllocatedVisAttributes = rhs.fAllocatedVisAttributes;
  rhs.fpVisAttributes = nullptr;
  rhs.fAllocatedVisAttributes = false;
  return *this;
}

void G4Visible::SetVisAttributes (const G4VisAttributes& VA) {
  // Allocate G4VisAttributes on the heap in case the user specifies a
  // short-lived VA for a long-lived G4Visible.  Flag so that it can
  // be deleted in the destructor.
  // First delete any G4VisAttributes already on the heap...
  if (fAllocatedVisAttributes) delete fpVisAttributes;
  fpVisAttributes = new G4VisAttributes(VA);
  fAllocatedVisAttributes = true;
}


void G4Visible::SetVisAttributes (const G4VisAttributes* pVA) {
  // First delete any G4VisAttributes already on the heap...
  if (fAllocatedVisAttributes) delete fpVisAttributes;
  fpVisAttributes = pVA;
  fAllocatedVisAttributes = false;
}

G4bool G4Visible::operator != (const G4Visible& right) const {
  if (fInfo != right.fInfo) return false;
  if ((fpVisAttributes != nullptr) && (right.fpVisAttributes != nullptr))
    return *fpVisAttributes != *right.fpVisAttributes;
  if ((fpVisAttributes == nullptr) && (right.fpVisAttributes == nullptr)) return false;
  return true;
}

std::ostream& operator << (std::ostream& os, const G4Visible& v) {
  os << "G4Visible: ";
  if (!v.fInfo.empty()) os << "User information: " << v.fInfo;
  os << '\n';
  if (v.fpVisAttributes != nullptr) return os << *(v.fpVisAttributes);
  return os << "No Visualization Attributes";
}
