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

#include "G4AttDef.hh"

#include "G4AttDefStore.hh"

using namespace std;

// Deprecated - see header.
std::ostream& operator<<
(std::ostream& os, const std::map<G4String,G4AttDef>* definitions)
{
  os << "G4AttDef: Deprecated output function.  Use const reference instead." << endl;
  if (definitions != nullptr) {
    os << *definitions;
  } else {
    os << "G4AttCheck: ERROR: zero definitions pointer." << endl;
  }
  return os;
}

std::ostream& operator<<
(std::ostream& os, const std::map<G4String,G4AttDef>& definitions)
{
  G4String storeKey;
  if (G4AttDefStore::GetStoreKey(&definitions, storeKey)) {
    os << storeKey << ":";
  }
  std::map<G4String,G4AttDef>::const_iterator i;
  for (i = definitions.begin(); i != definitions.end(); ++i) {
    const G4String& name = i->first;
    const G4AttDef& attDef = i->second;
    if (attDef.GetCategory() == "Physics") {
      os << "\n  " << attDef.GetDesc()
             << " (" << name << "): ";
      if (!attDef.GetExtra().empty()) {
        if (attDef.GetExtra() != "G4BestUnit") os << "unit: ";
        os << attDef.GetExtra() << " (";
      }
      os << attDef.GetValueType();
      if (!attDef.GetExtra().empty()) os << ")";
    }
  }
  os << endl;
  return os;
}
