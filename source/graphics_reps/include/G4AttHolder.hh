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
// $Id: G4AttHolder.hh,v 1.2 2006-10-24 06:00:29 allison Exp $
// GEANT4 tag $Name: not supported by cvs2svn $

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

class G4AttValue;
class G4AttDef;

class G4AttHolder {
public:
  G4AttHolder() {}
  /// Note: G4AttValues are deleted here...
  ~G4AttHolder();
  const std::vector<const std::vector<G4AttValue>*>& GetAttValues() const
  {return fValues;}
  const std::vector<const std::map<G4String,G4AttDef>*>& GetAttDefs() const
  {return fDefs;}
  const G4String& G4AttHolder::GetAttDefsName(size_t i) const;
  /// Add expendable G4AttValues - they will be deleted in destructor...
  void AddAttValues(const std::vector<G4AttValue>* values)
  {fValues.push_back(values);}
  /// G4AttDefs are assumed to have long life...
  void AddAttDefs(const std::map<G4String,G4AttDef>* defs)
  {fDefs.push_back(defs);}
private:
  std::vector<const std::vector<G4AttValue>*> fValues;
  std::vector<const std::map<G4String,G4AttDef>*> fDefs;
};

#endif //G4ATTHOLDER_HH
