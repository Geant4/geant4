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

#ifndef G4XNNstarTable_h
#define G4XNNstarTable_h

#include "globals.hh"
#include "G4PhysicsVector.hh"
#include "G4VXResonanceTable.hh"

#include <map>

class G4XNNstarTable 
{

public:

  // Constructor
  G4XNNstarTable(); 

  // Destructor
  virtual ~G4XNNstarTable();

  // Cross section table
  virtual const G4PhysicsVector* CrossSectionTable(const G4String& particleName) const;

  // Operators
  G4bool operator==(const G4XNNstarTable &right) const;
  G4bool operator!=(const G4XNNstarTable &right) const;


protected:


private:  

  G4XNNstarTable(const G4XNNstarTable &right);
  G4XNNstarTable& operator=(const G4XNNstarTable &right);

  std::map <G4String, G4double*, std::less<G4String> > xMap;

  static const G4int sizeNNstar;

  // The energies corresponding to the following cross sections
  static const G4double energyTable[121];

  // Cross sections for p p -> N N*
  static const G4double sigmaNN1440[121];
  static const G4double sigmaNN1520[121];
  static const G4double sigmaNN1535[121];
  static const G4double sigmaNN1650[121];
  static const G4double sigmaNN1675[121];
  static const G4double sigmaNN1680[121];
  static const G4double sigmaNN1700[121];
  static const G4double sigmaNN1710[121];
  static const G4double sigmaNN1720[121];
  static const G4double sigmaNN1900[121];
  static const G4double sigmaNN1990[121];
  static const G4double sigmaNN2090[121];
  static const G4double sigmaNN2190[121];
  static const G4double sigmaNN2220[121];
  static const G4double sigmaNN2250[121];

};

#endif

