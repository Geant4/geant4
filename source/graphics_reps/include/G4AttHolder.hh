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

#ifndef G4ATTHOLDER_HH
#define G4ATTHOLDER_HH

// Class Description:
//
// @class G4AttHolder
//
// @brief Holds G4AttValues and their corresponding G4AttDef map.
//
// Intended for inheritance by classes that need to hold them.
//
// For further details, see the HepRep home page at http://heprep.freehep.org
//  
// @author J.Allison
// @author J.Perl
// Class Description - End:

#include "globals.hh"

#include <vector>
#include <map>

#include "G4AttValue.hh"
#include "G4AttDef.hh"

class G4AttHolder {

public:

  G4AttHolder() = default;
  ~G4AttHolder();  // Note: G4AttValues are deleted here.
  G4AttHolder(const G4AttHolder&) = delete;  // Copy construction not allowed.
  G4AttHolder& operator=(const G4AttHolder&) = delete;  // Assignment not allowed.

  const std::vector<const std::vector<G4AttValue>*>& GetAttValues() const
  {return fValues;}
  const std::vector<const std::map<G4String,G4AttDef>*>& GetAttDefs() const
  {return fDefs;}

  void AddAtts(const std::vector<G4AttValue>* values,
               const std::map<G4String,G4AttDef>* defs)
  {fValues.push_back(values); fDefs.push_back(defs);}
  // Note: G4AttValues are assumed to be expendable - they will be
  // deleted in the destructor.  G4AttDefs are assumed to have long
  // life.

private:
  std::vector<const std::vector<G4AttValue>*> fValues;
  std::vector<const std::map<G4String,G4AttDef>*> fDefs;
};

#endif //G4ATTHOLDER_HH
