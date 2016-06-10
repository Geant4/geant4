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
// $Id: G3MatTable.hh 67982 2013-03-13 10:36:03Z gcosmo $
//
// ----------------------
// Class description:
//
// G3 materials table.
// Maps G3 material indices to their G4Material object counterparts.

// ----------------------
//
// by I.Hrivnacova, 27 Sep 99

#ifndef G3MATTABLE_HH
#define G3MATTABLE_HH 1

#include "G3MatTableEntry.hh"
#include "G3toG4Defs.hh"

#include "globals.hh"

#include <vector>

class G4Material;

typedef std::vector<G3MatTableEntry*>  G3MaterialVector;

class G3MatTable
{
  public: // with description

    G3MatTable();
    virtual ~G3MatTable();
    
    // methods
    G4Material* get(G4int id) const;
    void put(G4int id, G4Material* material);
    void Clear();

  private:

    G3MaterialVector*  fMatVector;
};

extern G3G4DLL_API G3MatTable G3Mat;

#endif
