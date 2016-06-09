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
// p p -> Delta Delta cross section tables
//
// -------------------------------------------------------------------

#ifndef G4XDeltaDeltaTable_h
#define G4XDeltaDeltaTable_h

#include "globals.hh"
#include "G4PhysicsVector.hh"
#include "G4PhysicsFreeVector.hh"
#include "G4VXResonanceTable.hh"


class G4XDeltaDeltaTable : public G4VXResonanceTable
{

public:

  G4XDeltaDeltaTable(); 

  ~G4XDeltaDeltaTable();

  virtual G4PhysicsVector* CrossSectionTable() const;

  G4bool operator==(const G4XDeltaDeltaTable &right) const;
  G4bool operator!=(const G4XDeltaDeltaTable &right) const;


protected:


private:  

  G4XDeltaDeltaTable(const G4XDeltaDeltaTable &right);
  G4XDeltaDeltaTable& operator=(const G4XDeltaDeltaTable &right);

  static const G4int sizeDeltaDelta;

  // The energies corresponding to the following cross sections
  static const G4double energyTable[121];

  // Cross sections for p p -> N N*
  static const G4double sigmaDD1232[121];

  G4int size;

};

#endif

