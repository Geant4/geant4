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

#ifndef G4BaryonWidth_h
#define G4BaryonWidth_h

#include "globals.hh"
#include "G4ResonanceWidth.hh"
#include <map>
#include <vector>

class G4PhysicsVector;

class G4BaryonWidth : public G4ResonanceWidth
{
public:

  G4BaryonWidth();

  virtual ~G4BaryonWidth();

  G4bool operator==(const G4BaryonWidth& right) const;
  G4bool operator!=(const G4BaryonWidth& right) const;

  // Returned pointer is owned by the client
  virtual G4PhysicsVector* MassDependentWidth(const G4String& name) const;

protected:

  
private:  

  G4BaryonWidth(const G4BaryonWidth& right);
  G4BaryonWidth& operator=(const G4BaryonWidth& right);

  // Map of PhysicsVectors for width interpolation
  std::map<G4String, G4double*, std::less<G4String> > wMap;
  //  std::map<G4String, G4double*, std::less<const G4String> > wMap;

  static const G4double baryonEnergyTable[120];

  static const G4double wN1440[120];
  static const G4double wN1520[120];
  static const G4double wN1535[120];
  static const G4double wN1650[120];
  static const G4double wN1675[120];
  static const G4double wN1680[120];
  static const G4double wN1700[120];
  static const G4double wN1710[120];
  static const G4double wN1720[120];
  static const G4double wN1900[120];
  static const G4double wN1990[120];
  static const G4double wN2090[120];
  static const G4double wN2190[120];
  static const G4double wN2220[120];
  static const G4double wN2250[120];

  static const G4double wD1600[120];
  static const G4double wD1620[120];
  static const G4double wD1700[120];
  static const G4double wD1900[120];
  static const G4double wD1905[120];
  static const G4double wD1910[120];
  static const G4double wD1920[120];
  static const G4double wD1930[120];
  static const G4double wD1950[120];

  static const G4double wDelta[120];

  static const G4double wL1405[120];
  static const G4double wL1520[120];
  static const G4double wL1600[120];
  static const G4double wL1670[120];
  static const G4double wL1690[120];
  static const G4double wL1800[120];
  static const G4double wL1810[120];
  static const G4double wL1820[120];
  static const G4double wL1830[120];
  static const G4double wL1890[120];
  static const G4double wL2100[120];
  static const G4double wL2110[120];

  static const G4double wS1385[120];
  static const G4double wS1660[120];
  static const G4double wS1670[120];
  static const G4double wS1750[120];
  static const G4double wS1775[120];
  static const G4double wS1915[120];
  static const G4double wS1940[120];
  static const G4double wS2030[120];

  static const G4double wX1530[120];
  static const G4double wX1690[120];
  static const G4double wX1820[120];
  static const G4double wX1950[120];
  static const G4double wX2030[120];

  G4int wSize;
};
  
#endif
