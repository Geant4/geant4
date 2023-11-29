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

#ifndef G4ATTCHECK_HH
#define G4ATTCHECK_HH

// Class Description:
//
// @class G4AttCheck
//
// @brief Checks G4AttValue's and their corresponding G4AttDef map.
//
// Usage (a): G4AttCheck(values,definitions);
//    or (b): G4cout << G4AttCheck(values,definitions) << G4endl;
//
// For further details, see the HepRep home page at http://heprep.freehep.org
//  
// @author J.Allison
// @author J.Perl
// Class Description - End:

#include "globals.hh"

#include <set>
#include <vector>
#include <map>
#include <iostream>

class G4AttValue;
class G4AttDef;

class G4AttCheck {

public: // With description

  G4AttCheck
  (const std::vector<G4AttValue>* values,
   const std::map<G4String,G4AttDef>* definitions);

  ~G4AttCheck() = default;

  const std::vector<G4AttValue>* GetAttValues() const {
    return fpValues;
  }

  const std::map<G4String,G4AttDef>* GetAttDefs() const {
    return fpDefinitions;
  }

  G4bool Check(const G4String& leader = "") const;
  // Check only.  Silent unless error - then G4cerr.  Returns error.
  // Provide leader, e.g., name of calling function.

  G4bool Standard
  (std::vector<G4AttValue>* standardValues,
   std::map<G4String,G4AttDef>* standardDefinitions)
    const;
  // Places standard versions in provided containers and returns error.

  friend std::ostream& operator<< (std::ostream&, const G4AttCheck&);

private:

  void AddValuesAndDefs
  (std::vector<G4AttValue>* newValues,
   std::map<G4String,G4AttDef>* newDefinitions,
   const G4String& oldName,
   const G4String& name,
   const G4String& value,
   const G4String& extra,
   const G4String& description = "") const;   // Utility function for Standard.

  void Init();   // Initialises maps and sets

  const std::vector<G4AttValue>* fpValues;
  const std::map<G4String,G4AttDef>* fpDefinitions;

  static G4ThreadLocal G4bool fFirst;  // Flag for initialising the following containers.
  static G4ThreadLocal std::set<G4String> *fUnitCategories;  // Set of legal unit categories.
  static G4ThreadLocal std::map<G4String,G4String> *fStandardUnits;  // Standard units.
  static G4ThreadLocal std::set<G4String> *fCategories;      // Set of legal categories.
  static G4ThreadLocal std::set<G4String> *fUnits;           // Set of legal units.
  static G4ThreadLocal std::set<G4String> *fValueTypes;      // Set of legal value types.
};

#endif //G4ATTCHECK_HH
