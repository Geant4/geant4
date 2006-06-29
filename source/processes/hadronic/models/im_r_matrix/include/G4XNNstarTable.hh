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

