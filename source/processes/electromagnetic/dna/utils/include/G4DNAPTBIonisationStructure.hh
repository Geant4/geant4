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
// Authors: S. Meylan and C. Villagrasa (IRSN, France)
// Models come from
// M. Bug et al, Rad. Phys and Chem. 130, 459-479 (2017)

#ifndef G4DNAPTBIonisationStructure_HH
#define G4DNAPTBIonisationStructure_HH 1
#include "globals.hh"
#include <map>
#include <vector>
class G4Material;
class G4DNAPTBIonisationStructure
{
 public:
  G4DNAPTBIonisationStructure();
  ~G4DNAPTBIonisationStructure() = default;
  G4double IonisationEnergy(G4int level, const size_t& materialName);
  G4int NumberOfLevels(const size_t& materialName);
  G4DNAPTBIonisationStructure(const G4DNAPTBIonisationStructure&) = delete;  // prevent copy-construction
  G4DNAPTBIonisationStructure& operator=(
    const G4DNAPTBIonisationStructure& right) = delete;  // prevent assignement

 private:
  // Number of Ionisation levels of the water molecule
  std::map<size_t, G4int> nLevels;
  std::map<size_t, std::vector<G4double>> energyConstant;
  size_t ReplaceMaterial(const size_t& materialID);
  G4Material* fpGuanine_PU = nullptr;
  G4Material* fpTHF = nullptr;
  G4Material* fpPY = nullptr;
  G4Material* fpPU = nullptr;
  G4Material* fpTMP = nullptr;
  G4Material* fpG4_WATER = nullptr;
  G4Material* fpBackbone_THF = nullptr;
  G4Material* fpCytosine_PY = nullptr;
  G4Material* fpThymine_PY = nullptr;
  G4Material* fpAdenine_PU = nullptr;
  G4Material* fpBackbone_TMP = nullptr;
  G4Material* fpN2 = nullptr;
};

#endif
