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
// p p -> N Delta* cross section tables
//
// -------------------------------------------------------------------

#ifndef G4XDeltaDeltastarTable_h
#define G4XDeltaDeltastarTable_h

#include "globals.hh"
#include "G4PhysicsVector.hh"
#include "G4VXResonanceTable.hh"

#include <map>

class G4XDeltaDeltastarTable 
{

public:

  G4XDeltaDeltastarTable(); 

  virtual ~G4XDeltaDeltastarTable();

  virtual const G4PhysicsVector* CrossSectionTable(const G4String& particleName) const;

  G4bool operator==(const G4XDeltaDeltastarTable &right) const;
  G4bool operator!=(const G4XDeltaDeltastarTable &right) const;


protected:


private:  

  G4XDeltaDeltastarTable(const G4XDeltaDeltastarTable &right);
  G4XDeltaDeltastarTable& operator=(const G4XDeltaDeltastarTable &right);

  std::map <G4String, G4double*, std::less<G4String> > xMap;

  static const G4int sizeDeltaDeltastar;

  // The energies corresponding to the following cross sections
  static const G4double energyTable[121];

  // Cross sections for p p -> N Delta*
  static const G4double sigmaDD1600[121];
  static const G4double sigmaDD1620[121];
  static const G4double sigmaDD1700[121];
  static const G4double sigmaDD1900[121];
  static const G4double sigmaDD1905[121];
  static const G4double sigmaDD1910[121];
  static const G4double sigmaDD1920[121];
  static const G4double sigmaDD1930[121]; // 40 is missing... @@@@@@@
  static const G4double sigmaDD1950[121];

};

#endif

