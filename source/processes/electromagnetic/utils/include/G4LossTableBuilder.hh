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
//
// -------------------------------------------------------------------
//
// GEANT4 Class header file
//
//
// File name:     G4LossTableBuilder
//
// Author:        Vladimir Ivanchenko on base of Laszlo Urban code
// 
// Creation date: 03.01.2002
//
// Modifications: 
//
//
// Class Description: 
//
// Provide building of dE/dx, range, and inverse range tables.

// -------------------------------------------------------------------
//

#ifndef G4LossTableBuilder_h
#define G4LossTableBuilder_h 1

#include "g4std/vector"
#include "globals.hh"
#include "G4PhysicsTable.hh"


class G4LossTableBuilder 
{

public:

  G4LossTableBuilder() {};

  ~G4LossTableBuilder() {};

  G4PhysicsTable* BuildDEDXTable(const G4std::vector<G4PhysicsTable*>&);

  G4PhysicsTable* BuildRangeTable(const G4PhysicsTable* dedxTable);

  G4PhysicsTable* BuildInverseRangeTable(const G4PhysicsTable* dedxTable,
                                         const G4PhysicsTable* rangeTable);
 
private:

  G4LossTableBuilder & operator=(const  G4LossTableBuilder &right);
  G4LossTableBuilder(const  G4LossTableBuilder&);

};

//....oooOO0OOooo.......oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
