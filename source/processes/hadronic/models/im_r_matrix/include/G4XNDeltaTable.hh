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
// Hadron Kinetic Model
// p p -> Delta Delta cross section tables
//
// -------------------------------------------------------------------

#ifndef G4XNDeltaTable_h
#define G4XNDeltaTable_h

#include "globals.hh"
#include "G4PhysicsVector.hh"
#include "G4PhysicsFreeVector.hh"
#include "G4VXResonanceTable.hh"


class G4XNDeltaTable : public G4VXResonanceTable
{

public:

  G4XNDeltaTable(); 

  ~G4XNDeltaTable();

  virtual G4PhysicsVector* CrossSectionTable() const;

  G4bool operator==(const G4XNDeltaTable &right) const;
  G4bool operator!=(const G4XNDeltaTable &right) const;


protected:


private:  

  G4XNDeltaTable(const G4XNDeltaTable &right);
  G4XNDeltaTable& operator=(const G4XNDeltaTable &right);

  // The energies corresponding to the following cross sections
  static const G4double energyTable[121];

  // Cross sections for p p -> N N*
  static const G4double sigmaND1232[121];

  G4int size;

};

#endif

