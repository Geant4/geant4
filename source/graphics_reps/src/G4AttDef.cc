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
// $Id$

#include "G4AttDef.hh"

#include "G4AttDefStore.hh"

std::ostream& operator<<
  (std::ostream& os, const std::map<G4String,G4AttDef>* definitions)
{
  using namespace std;
  if (!definitions) {
    os << "G4AttCheck: ERROR: zero definitions pointer." << endl;
    return os;
  }
  G4String storeKey;
  if (G4AttDefStore::GetStoreKey(definitions, storeKey)) {
    os << storeKey << ":";
  }
  std::map<G4String,G4AttDef>::const_iterator i;
  for (i = definitions->begin(); i != definitions->end(); ++i) {
    if (i->second.GetCategory() == "Physics") {
      os << "\n  " << i->second.GetDesc()
	     << " (" << i->first << "): ";
      if (!i->second.GetExtra().empty()) {
	if (i->second.GetExtra() != "G4BestUnit") os << "unit: ";
	os << i->second.GetExtra() << " (";
      }
      os << i->second.GetValueType();
      if (!i->second.GetExtra().empty()) os << ")";
    }
  }
  os << endl;
  return os;
}
