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
// -------------------------------------------------------------------
//      GEANT4 Class file
//
//
//      File name:     G4BaryonPartialWidth
//
//      Author:        Maria Grazia Pia (MariaGrazia.Pia@genova.infn.it)
// 
//      Creation date: 15 April 1999
//
//      Modifications: 
//      
// -------------------------------------------------------------------

#ifndef G4BARYONPARTIALWIDTH_HH
#define G4BARYONPARTIALWIDTH_HH

#include "globals.hh"
#include "G4ResonancePartialWidth.hh"
#include <map>

class G4PhysicsVector;

class G4BaryonPartialWidth :public G4ResonancePartialWidth
{
public:

  G4BaryonPartialWidth();

  virtual ~G4BaryonPartialWidth();

  // Returned pointer is owned by the client
  virtual G4PhysicsVector* MassDependentWidth(const G4String& name) const;

protected:
  
private:  

  G4BaryonPartialWidth(const G4BaryonPartialWidth& right);
  G4BaryonPartialWidth& operator=(const G4BaryonPartialWidth& right);

  std::map<G4String, G4double*, std::less<G4String> > wMap;

  static const G4double energies[120];

  static const G4double pwN1440_Npi[120];
  static const G4double pwN1440_Npipi[120];
  static const G4double pwN1440_Dpi[120];

  static const G4double pwN1520_Ngamma[120];
  static const G4double pwN1520_Npi[120];
  static const G4double pwN1520_Npipi[120];
  static const G4double pwN1520_Dpi[120];

  static const G4double pwN1535_Ngamma[120];
  static const G4double pwN1535_Npi[120];
  static const G4double pwN1535_Neta[120];
  static const G4double pwN1535_Npipi[120];
  static const G4double pwN1535_Nstarpi[120];

  static const G4double pwN1650_Ngamma[120];
  static const G4double pwN1650_Npi[120];
  static const G4double pwN1650_Neta[120] ;
  static const G4double pwN1650_Npipi[120];
  static const G4double pwN1650_Dpi[120];
  static const G4double pwN1650_Nstarpi[120];
  static const G4double pwN1650_LK[120];

  static const G4double pwN1675_Npi[120];
  static const G4double pwN1675_Dpi[120];

  static const G4double pwN1680_Ngamma[120];
  static const G4double pwN1680_Npi[120];
  static const G4double pwN1680_Npipi[120];
  static const G4double pwN1680_Dpi[120];
 
  static const G4double pwN1700_Npi[120];
  static const G4double pwN1700_Neta[120];
  static const G4double pwN1700_Nrho[120];
  static const G4double pwN1700_Npipi[120];
  static const G4double pwN1700_Dpi[120];

  static const G4double pwN1710_Npi[120];
  static const G4double pwN1710_Neta[120];
  static const G4double pwN1710_Nrho[120];
  static const G4double pwN1710_Npipi[120];
  static const G4double pwN1710_Dpi[120];
  static const G4double pwN1710_Nstarpi[120];
  static const G4double pwN1710_LK[120];

  static const G4double pwN1720_Ngamma[120];
  static const G4double pwN1720_Npi[120];
  static const G4double pwN1720_Nrho[120];
  static const G4double pwN1720_Npipi[120];
  static const G4double pwN1720_Dpi[120];
  static const G4double pwN1720_LK[120];

  static const G4double pwN1900_Npi[120];
  static const G4double pwN1900_Nomega[120];
  static const G4double pwN1900_Nrho[120];
  static const G4double pwN1900_Dpi[120];

  static const G4double pwN1990_Npi[120];
  static const G4double pwN1990_Nrho[120];
  static const G4double pwN1990_Npipi[120];
  static const G4double pwN1990_Dpi[120];
  static const G4double pwN1990_Nstarpi[120];
  static const G4double pwN1990_LK[120];

  static const G4double pwN2090_Npi[120];
  static const G4double pwN2090_Neta[120];
  static const G4double pwN2090_Nrho[120];
  static const G4double pwN2090_Npipi[120];
  static const G4double pwN2090_Dpi[120];

  static const G4double pwN2190_Npi[120];
  static const G4double pwN2190_Nrho[120];
  static const G4double pwN2190_Npipi[120];
  static const G4double pwN2190_Dpi[120];
  static const G4double pwN2190_Nstarpi[120];

  static const G4double pwN2220_Npi[120];
  static const G4double pwN2220_Nrho[120];
  static const G4double pwN2220_Npipi[120];
  static const G4double pwN2220_Dpi[120];
  static const G4double pwN2250_Npi[120];
  static const G4double pwN2250_Nrho[120];
  static const G4double pwN2250_Npipi[120];
  static const G4double pwN2250_Dpi[120];
  static const G4double pwN2250_Nstarpi[120];

  static const G4double pwD1232_Ngamma[120];
  static const G4double pwD1232_Npi[120];

  static const G4double pwD1600_Npi[120];
  static const G4double pwD1600_Dpi[120];
  static const G4double pwD1600_Nstarpi[120];

  static const G4double pwD1620_Ngamma[120];
  static const G4double pwD1620_Npi[120];
  static const G4double pwD1620_Dpi[120];
  static const G4double pwD1620_Nstarpi[120];

  static const G4double pwD1700_Ngamma[120];
  static const G4double pwD1700_Npi[120];
  static const G4double pwD1700_Nrho[120];
  static const G4double pwD1700_Dpi[120];
  static const G4double pwD1700_Nstarpi[120];

  static const G4double pwD1900_Npi[120];
  static const G4double pwD1900_Nrho[120];
  static const G4double pwD1900_Dpi[120];
  static const G4double pwD1900_Nstarpi[120];

  static const G4double pwD1905_Ngamma[120];
  static const G4double pwD1905_Npi[120];
  static const G4double pwD1905_Nrho[120];
  static const G4double pwD1905_Dpi[120];
  static const G4double pwD1905_Nstarpi[120];

  static const G4double pwD1910_Npi[120];
  static const G4double pwD1910_Nrho[120];
  static const G4double pwD1910_Dpi[120];
  static const G4double pwD1910_Nstarpi[120];

  static const G4double pwD1920_Npi[120];
  static const G4double pwD1920_Nrho[120];
  static const G4double pwD1920_Dpi[120];
  static const G4double pwD1920_Nstarpi[120];

  static const G4double pwD1930_Npi[120];
  static const G4double pwD1930_Nrho[120];
  static const G4double pwD1930_Dpi[120];
  static const G4double pwD1930_Nstarpi[120];

  static const G4double pwD1950_Ngamma[120];
  static const G4double pwD1950_Npi[120];
  static const G4double pwD1950_Nrho[120];
  static const G4double pwD1950_Dpi[120];
  static const G4double pwD1950_Nstarpi[120];

  static const G4double pwL1405_Spi[120];

  static const G4double pwL1520_NKbar[120];
  static const G4double pwL1520_Spi[120];
  static const G4double pwL1520_Sstarpi[120];
  static const G4double pwL1520_Lgamma[120];

  static const G4double pwL1600_NKbar[120];
  static const G4double pwL1600_Spi[120];
  static const G4double pwL1670_NKbar[120];
  static const G4double pwL1670_Spi[120];
  static const G4double pwL1670_Leta[120];

  static const G4double pwL1690_NKbar[120];
  static const G4double pwL1690_Spi[120];
  static const G4double pwL1690_Sstarpi[120];
  
  static const G4double pwL1800_NKbar[120];
  static const G4double pwL1800_NKstarbar[120];
  static const G4double pwL1800_Spi[120];
  static const G4double pwL1800_Sstarpi[120];

  static const G4double pwL1810_NKbar[120];
  static const G4double pwL1810_NKstarbar[120];
  static const G4double pwL1810_Spi[120];
  static const G4double pwL1810_Sstarpi[120];

  static const G4double pwL1820_NKbar[120];
  static const G4double pwL1820_Spi[120];
  static const G4double pwL1820_Sstarpi[120];

  static const G4double pwL1830_NKbar[120];
  static const G4double pwL1830_Spi[120];
  static const G4double pwL1830_Sstarpi[120];

  static const G4double pwL1890_NKbar[120];
  static const G4double pwL1890_NKstarbar[120];
  static const G4double pwL1890_Spi[120];
  static const G4double pwL1890_Sstarpi[120];

  static const G4double pwL2100_NKbar[120];
  static const G4double pwL2100_NKstarbar[120];
  static const G4double pwL2100_Spi[120];
  static const G4double pwL2100_Sstarpi[120];
  static const G4double pwL2100_Leta[120];
  static const G4double pwL2100_Lomega[120];
  static const G4double pwL2110_NKbar[120];
  static const G4double pwL2110_NKstarbar[120];
  static const G4double pwL2110_Spi[120];

  static const G4double pwS1385_Spi[120];
  static const G4double pwS1385_Lpi[120];

  static const G4double pwS1660_NKbar[120];
  static const G4double pwS1660_Spi[120];
  static const G4double pwS1660_Lpi[120];

  static const G4double pwS1670_NKbar[120];
  static const G4double pwS1670_Spi[120];
  static const G4double pwS1670_Lpi[120];

  static const G4double pwS1750_NKbar[120];
  static const G4double pwS1750_Spi[120];
  static const G4double pwS1750_Seta[120];

  static const G4double pwS1775_NKbar[120];
  static const G4double pwS1775_Spi[120];
  static const G4double pwS1775_Sstarpi[120];
  static const G4double pwS1775_Lpi[120];
  static const G4double pwS1775_Lstarpi[120];

  static const G4double pwS1915_NKbar[120];
  static const G4double pwS1915_Spi[120];
  static const G4double pwS1915_Sstarpi[120];
  static const G4double pwS1915_Lpi[120];

  static const G4double pwS1940_NKbar[120];
  static const G4double pwS1940_NKstarbar[120];
  static const G4double pwS1940_Spi[120];
  static const G4double pwS1940_Sstarpi[120];
  static const G4double pwS1940_Lpi[120];
  static const G4double pwS1940_Lstarpi[120];
  static const G4double pwS1940_DKbar[120];

  static const G4double pwS2030_NKbar[120];
  static const G4double pwS2030_NKstarbar[120];
  static const G4double pwS2030_Spi[120];
  static const G4double pwS2030_Sstarpi[120];
  static const G4double pwS2030_Lpi[120];
  static const G4double pwS2030_Lstarpi[120];
  static const G4double pwS2030_DKbar[120];

  static const G4double pwX1530_Xpi[120];
  static const G4double pwX1530_Xgamma[120];

  static const G4double pwX1690_Xpi[120];
  static const G4double pwX1690_LKbar[120];
  static const G4double pwX1690_SKbar[120];

  static const G4double pwX1820_Xpi[120];
  static const G4double pwX1820_LKbar[120];
  static const G4double pwX1820_SKbar[120];

  static const G4double pwX1950_Xpi[120];
  static const G4double pwX1950_LKbar[120];
  static const G4double pwX1950_SKbar[120];

  static const G4double pwX2030_Xpi[120];
  static const G4double pwX2030_LKbar[120];
  static const G4double pwX2030_SKbar[120];

  G4int wSize;

};
  
#endif
