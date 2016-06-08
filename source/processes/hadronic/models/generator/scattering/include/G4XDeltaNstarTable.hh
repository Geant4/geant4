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

#ifndef G4XDeltaNstarTable_h
#define G4XDeltaNstarTable_h

#include "globals.hh"
#include "G4PhysicsVector.hh"
#include "G4VXResonanceTable.hh"

#include "g4std/map"

class G4XDeltaNstarTable 
{

public:

  G4XDeltaNstarTable(); 

  virtual ~G4XDeltaNstarTable();

  virtual const G4PhysicsVector* CrossSectionTable(const G4String& particleName) const;

  G4bool operator==(const G4XDeltaNstarTable &right) const;
  G4bool operator!=(const G4XDeltaNstarTable &right) const;


protected:

private:  

  G4XDeltaNstarTable(const G4XDeltaNstarTable &right);
  G4XDeltaNstarTable& operator=(const G4XDeltaNstarTable &right);

  G4std::map <G4String, G4double*, G4std::less<G4String> > xMap;
  
  static const G4int sizeDeltaNstar;

  // The energies corresponding to the following cross sections
  static const G4double energyTable[121];

  // Cross sections for p p -> N N*
  static const G4double sigmaDN1440[121];
  static const G4double sigmaDN1520[121];
  static const G4double sigmaDN1535[121];
  static const G4double sigmaDN1650[121];
  static const G4double sigmaDN1675[121];
  static const G4double sigmaDN1680[121];
  static const G4double sigmaDN1700[121];
  static const G4double sigmaDN1710[121];
  static const G4double sigmaDN1720[121];
  static const G4double sigmaDN1900[121];
  static const G4double sigmaDN1990[121];
  static const G4double sigmaDN2090[121];
  static const G4double sigmaDN2190[121];
  static const G4double sigmaDN2220[121];
  static const G4double sigmaDN2250[121];

};

#endif

