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
// G4tgbMaterialPropertiesTable
//
// Class description:
//
// Transient class of a material property table; 
// builds a G4tgbMaterialPropertiesTable.
//
#ifndef G4tgbMaterialPropertiesTable_hh
#define G4tgbMaterialPropertiesTable_hh 1

#include <vector>
#include <string>

#include "globals.hh"
#include "G4tgrMaterialPropertiesTable.hh"
#include "G4MaterialPropertiesTable.hh"

class G4tgbMaterialPropertiesTable
{
  public:

    G4tgbMaterialPropertiesTable(G4tgrMaterialPropertiesTable*);
    ~G4tgbMaterialPropertiesTable();

    G4String GetName() const { return theName; }
    G4MaterialPropertiesTable* BuildG4MaterialPropertiesTable();
      // Build a G4MaterialPropertiesTable

  private:

    G4String theName;
    G4tgrMaterialPropertiesTable* theTgrTable = nullptr;
    G4MaterialPropertiesTable* theTable = new G4MaterialPropertiesTable();
};

#endif

