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
// $Id: G3MedTable.hh 67982 2013-03-13 10:36:03Z gcosmo $
//
// ----------------------
// Class description:
//
// G3 tracking medium table.
// Maps G3 tracking medium indices to their corresponding
// G4 objects (G4Material, G4MagneticField, G4UserLimits)
// and the sensitivity flag (isvol).

// ----------------------
//
// by I.Hrivnacova, 27 Sep 99

#ifndef G3MEDTABLE_HH
#define G3MEDTABLE_HH 1

#include "G3MedTableEntry.hh"
#include "G3toG4Defs.hh"

#include "globals.hh"

#include <vector>

class G4Material;
class G4MagneticField;
class G4UserLimits;

typedef std::vector<G3MedTableEntry*>  G3MediumVector;

class G3MedTable
{
  public: // with description

    G3MedTable();
    virtual ~G3MedTable();
    
    // methods
    G3MedTableEntry* get(G4int id) const;
    void put(G4int id, G4Material* material, G4MagneticField* field,
             G4UserLimits* limits, G4int isvol);
    G4int GetSize() const;
    G3MedTableEntry*  GetMTE(G4int i) const;
    void Clear();

  private:

    G3MediumVector*  fMedVector;
};

extern G3G4DLL_API G3MedTable G3Med;
#endif
