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

  static G4double baryonEnergyTable[120];

  static G4double wN1440[120];
  static G4double wN1520[120];
  static G4double wN1535[120];
  static G4double wN1650[120];
  static G4double wN1675[120];
  static G4double wN1680[120];
  static G4double wN1700[120];
  static G4double wN1710[120];
  static G4double wN1720[120];
  static G4double wN1900[120];
  static G4double wN1990[120];
  static G4double wN2090[120];
  static G4double wN2190[120];
  static G4double wN2220[120];
  static G4double wN2250[120];

  static G4double wD1600[120];
  static G4double wD1620[120];
  static G4double wD1700[120];
  static G4double wD1900[120];
  static G4double wD1905[120];
  static G4double wD1910[120];
  static G4double wD1920[120];
  static G4double wD1930[120];
  static G4double wD1950[120];

  static G4double wDelta[120];

  static G4double wL1405[120];
  static G4double wL1520[120];
  static G4double wL1600[120];
  static G4double wL1670[120];
  static G4double wL1690[120];
  static G4double wL1800[120];
  static G4double wL1810[120];
  static G4double wL1820[120];
  static G4double wL1830[120];
  static G4double wL1890[120];
  static G4double wL2100[120];
  static G4double wL2110[120];

  static G4double wS1385[120];
  static G4double wS1660[120];
  static G4double wS1670[120];
  static G4double wS1750[120];
  static G4double wS1775[120];
  static G4double wS1915[120];
  static G4double wS1940[120];
  static G4double wS2030[120];

  static G4double wX1530[120];
  static G4double wX1690[120];
  static G4double wX1820[120];
  static G4double wX1950[120];
  static G4double wX2030[120];

  G4int wSize;
};
  
#endif
