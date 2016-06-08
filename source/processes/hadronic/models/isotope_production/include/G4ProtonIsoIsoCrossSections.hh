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
#ifndef G4ProtonIsoIsoCrossSections_h
#define G4ProtonIsoIsoCrossSections_h

#include "G4NeutronHPVector.hh"
#include "G4IsoProdCrossSections.hh"
#include "G4NeutronHPNames.hh"

class G4ProtonIsoIsoCrossSections
{

public:

  G4ProtonIsoIsoCrossSections();
  ~G4ProtonIsoIsoCrossSections();
  void Init(G4int A, G4int Z, G4double frac);
  G4double GetCrossSection(G4double anEnergy);
  G4String GetProductIsotope(G4double anEnergy);
  
  G4int GetZ() { return theZ; }
  G4int GetA() { return theA; }
  
private:
  
  G4int theNumberOfProducts;
  G4IsoProdCrossSections ** theProductionData;
  G4NeutronHPVector theCrossSection;
  G4NeutronHPNames theNames;
  G4bool hasData;
  G4int theA;
  G4int theZ;
};

#endif
