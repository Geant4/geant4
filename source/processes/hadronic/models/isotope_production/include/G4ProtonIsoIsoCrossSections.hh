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
