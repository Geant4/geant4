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

#ifndef G4XDeltaNstarTable_h
#define G4XDeltaNstarTable_h

#include "globals.hh"
#include "G4PhysicsVector.hh"
#include "G4VXResonanceTable.hh"

#include <map>

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

  std::map <G4String, G4double*, std::less<G4String> > xMap;
  
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

