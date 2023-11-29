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

#ifndef G4DNASCAVENGERMATERIAL_HH
#define G4DNASCAVENGERMATERIAL_HH
#include "globals.hh"
#include "G4ios.hh"
#include <map>
#include <vector>
#include "G4MoleculeCounter.hh"
#include "G4VScavengerMaterial.hh"
class G4Material;
class G4MolecularConfiguration;
class G4VChemistryWorld;

class G4DNAScavengerMaterial : public G4VScavengerMaterial
{
 public:
  using NbMoleculeInTime =
    std::map<G4double, int64_t, G4::MoleculeCounter::TimePrecision>;
  using MolType            = const G4MolecularConfiguration*;
  using MaterialMap        = std::map<MolType, int64_t>;
  using ReactantList       = std::vector<MolType>;
  using CounterMapType     = std::map<MolType, NbMoleculeInTime>;
  G4DNAScavengerMaterial() = default;
  explicit G4DNAScavengerMaterial(G4VChemistryWorld*);
  ~G4DNAScavengerMaterial() override                          = default;
  G4DNAScavengerMaterial(const G4DNAScavengerMaterial& right) = delete;
  G4DNAScavengerMaterial& operator=(const G4DNAScavengerMaterial&) = delete;
  void Initialize();

  void ReduceNumberMoleculePerVolumeUnitForMaterialConf(MolType, G4double);
  void AddNumberMoleculePerVolumeUnitForMaterialConf(MolType, G4double);
  G4double GetNumberMoleculePerVolumeUnitForMaterialConf(MolType) const;

  void AddAMoleculeAtTime(MolType, G4double time,
                          const G4ThreeVector* position = nullptr,
                          G4int number                  = 1);
  void RemoveAMoleculeAtTime(MolType, G4double time,
                             const G4ThreeVector* position = nullptr,
                             G4int number                  = 1);

  void Reset() override;

  void PrintInfo();

  MaterialMap::iterator end() { return fScavengerTable.end(); }
  MaterialMap::iterator begin() { return fScavengerTable.begin(); }
  size_t size() { return fScavengerTable.size(); }

  G4bool find(MolType type)
  {
    auto it = fScavengerTable.find(type);
    if(it != fScavengerTable.end())
    {
      return it->second > 0;
    }
    else
    {
      return false;
    }
  }

  void SetCounterAgainstTime() { fCounterAgainstTime = true; }

  std::vector<MolType> GetScavengerList() const
  {
    std::vector<MolType> output;
    for(const auto& it : fScavengerTable)
    {
      output.push_back(it.first);
    }
    return output;
  }

  void Dump();
  int64_t GetNMoleculesAtTime(MolType molecule, G4double time);
  G4bool SearchTimeMap(MolType molecule);
  int64_t SearchUpperBoundTime(G4double time, G4bool sameTypeOfMolecule);

 private:
  G4VChemistryWorld* fpChemistryInfo;
  G4bool fIsInitialized;
  MaterialMap fScavengerTable;
  CounterMapType fCounterMap;
  G4bool fCounterAgainstTime;
  G4int fVerbose;
  struct Search
  {
    Search() { fLowerBoundSet = false; }
    CounterMapType::iterator fLastMoleculeSearched;
    NbMoleculeInTime::iterator fLowerBoundTime;
    G4bool fLowerBoundSet;
  };

  std::unique_ptr<Search> fpLastSearch;
};
#endif  // G4DNASCAVENGERMATERIAL_HH
