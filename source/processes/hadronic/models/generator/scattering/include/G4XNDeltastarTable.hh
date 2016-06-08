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
// p p -> N Delta* cross section tables
//
// -------------------------------------------------------------------

#ifndef G4XNDELTASTARTABLE_HH
#define G4XNDELTASTARTABLE_HH

#include "globals.hh"
#include "G4PhysicsVector.hh"
#include "G4VXResonanceTable.hh"

#include "g4std/map"

class G4XNDeltastarTable 
{

public:

  // Constructor
  G4XNDeltastarTable(); 

  // Destructor
  virtual ~G4XNDeltastarTable();

  // Cross section table
  virtual const G4PhysicsVector* CrossSectionTable(const G4String& particleName) const;

  G4bool operator==(const G4XNDeltastarTable &right) const;
  G4bool operator!=(const G4XNDeltastarTable &right) const;


protected:


private:  

  G4XNDeltastarTable(const G4XNDeltastarTable &right);
  G4XNDeltastarTable& operator=(const G4XNDeltastarTable &right);

  G4std::map <G4String, G4double*, G4std::less<G4String> > xMap;

  static const G4int sizeNDeltastar;

  // The energies corresponding to the following cross sections
  static const G4double energyTable[121];

  // Cross sections for p p -> N Delta*
  static const G4double sigmaND1600[121];
  static const G4double sigmaND1620[121];
  static const G4double sigmaND1700[121];
  static const G4double sigmaND1900[121];
  static const G4double sigmaND1905[121];
  static const G4double sigmaND1910[121];
  static const G4double sigmaND1920[121]; 
  static const G4double sigmaND1930[121]; // 40 is missing... @@@@@@@
  static const G4double sigmaND1950[121];

};

#endif

