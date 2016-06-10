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

  virtual ~G4XNDeltaTable();

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

