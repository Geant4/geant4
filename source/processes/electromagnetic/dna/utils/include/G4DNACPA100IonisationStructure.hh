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
// Based on the work of M. Terrissol and M. C. Bordage
//
// Users are requested to cite the following papers:
// - M. Terrissol, A. Baudre, Radiat. Prot. Dosim. 31 (1990) 175-177
// - M.C. Bordage, J. Bordes, S. Edel, M. Terrissol, X. Franceries,
//   M. Bardies, N. Lampe, S. Incerti, Phys. Med. 32 (2016) 1833-1840
//
// Authors of this class:
// M.C. Bordage, M. Terrissol, S. Edel, J. Bordes, S. Incerti
//
// 15.01.2014: creation
//
//
// Based on the study by S. Zein et. al. Nucl. Inst. Meth. B 488 (2021) 70-82
// 1/2/2023 : Hoang added modification for DNA cross sections

#ifndef G4DNACPA100IonisationStructure_hh
#define G4DNACPA100IonisationStructure_hh

#include "globals.hh"

#include <vector>

class G4Material;

class G4DNACPA100IonisationStructure
{
  public:
    G4DNACPA100IonisationStructure();

    ~G4DNACPA100IonisationStructure() = default;

    G4double IonisationEnergy(const std::size_t& level, const std::size_t& MatID);

    G4double UEnergy(const std::size_t& level, const std::size_t& MatID);

    inline std::size_t NumberOfLevels(const std::size_t& MatID)
    {
      return fnLevels[MatID];
    }

  private:
    void InitialiseGuanine();
    void InitialiseWater();
    void InitialiseDeoxyribose();
    void InitialiseCytosine();
    void InitialiseThymine();
    void InitialiseAdenine();
    void InitialisePhosphate();

    std::map<std::size_t /*MatIndex*/, std::size_t> fnLevels;
    std::map<std::size_t /*MatIndex*/, std::vector<G4double>> fEnergyConstant;
    std::map<std::size_t /*MatIndex*/, std::vector<G4double>> fUConstant;

    G4Material* fpGuanine = nullptr;
    G4Material* fpG4_WATER = nullptr;
    G4Material* fpDeoxyribose = nullptr;
    G4Material* fpCytosine = nullptr;
    G4Material* fpThymine = nullptr;
    G4Material* fpAdenine = nullptr;
    G4Material* fpPhosphate = nullptr;
};

#endif  // G4DNACPA100IonisationStructure_hh
