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
// $Id: G3MedTable.hh,v 1.11 2003-01-27 10:45:40 gcosmo Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
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

#include "globals.hh"

#include "g4std/vector"

class G4Material;
class G4MagneticField;
class G4UserLimits;

typedef G4std::vector<G3MedTableEntry*>  G3MediumVector;

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

extern G3MedTable G3Med;
#endif
