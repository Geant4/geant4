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

#ifndef G4DNAPTBExcitationStructure_h
#define G4DNAPTBExcitationStructure_h 1
 
#include "globals.hh"
#include <vector>
#include <map>
class G4Material;
class G4DNAPTBExcitationStructure
{
public:

  G4DNAPTBExcitationStructure();
  ~G4DNAPTBExcitationStructure() = default;
  G4double ExcitationEnergy(const G4int& ExcLevel, const size_t &materialID);
  G4int NumberOfExcLevels(const size_t& materialID);
  G4DNAPTBExcitationStructure(const G4DNAPTBExcitationStructure&) = delete;  // prevent copy-construction
  G4DNAPTBExcitationStructure& operator=(
    const G4DNAPTBExcitationStructure& right) = delete;  // prevent assignement
private:
  // Number of Excitation levels of the water molecule
  std::map<size_t, G4int> nExcLevels;
  std::map<size_t, std::vector<G4double> > energyConstant;
  size_t ReplaceMaterial(const size_t &materialName);
  G4Material* fpN2;
};

#endif
