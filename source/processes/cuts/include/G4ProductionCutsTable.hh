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
// G4ProductionCutsTable
//
// Class description:
//
// G4ProductionCutsTable is a singleton class for a table of
// G4ProductionCuts objects. This class manages tables of production
// cuts and energy cuts for each particle type.

// Author: M.Asai, 5 October 2002 - First implementation
// Modifications: H.Kurashige, 2004-2008
// --------------------------------------------------------------------
#ifndef G4ProductionCutsTable_hh 
#define G4ProductionCutsTable_hh 1

#include <cmath>
#include <vector>

#include "globals.hh"
#include "G4ios.hh"
#include "G4MaterialCutsCouple.hh"
#include "G4MCCIndexConversionTable.hh"
#include "G4Region.hh"

class G4RegionStore;
class G4VRangeToEnergyConverter;
class G4LogicalVolume;
class G4VPhysicalVolume;
class G4ProductionCuts;
class G4ProductionCutsTableMessenger;

class G4ProductionCutsTable  
{
  public:

    static G4ProductionCutsTable* GetProductionCutsTable();
      // This static method returns the singleton pointer of this class object.
      // At first invocation, the singleton object is instantiated

    G4ProductionCutsTable(const G4ProductionCutsTable&) = delete;
    G4ProductionCutsTable& operator=(const G4ProductionCutsTable&) = delete;

    virtual ~G4ProductionCutsTable();

    void CreateCoupleTables();
      // Creates material cuts couples table and allocate the other tables

    void UpdateCoupleTable(G4VPhysicalVolume* currentWorld);
      // Triggers an update of the table of G4ProductionCuts objects

    void SetEnergyRange(G4double lowedge, G4double highedge);
      // Sets the limits of energy cuts for all particles

    G4double GetLowEdgeEnergy() const;
    G4double GetHighEdgeEnergy() const;
      // Get the limits of energy cuts for all particles

    G4double GetMaxEnergyCut();
    void SetMaxEnergyCut(G4double value);
      // Get/set max cut energy of RangeToEnergy converter
      // for all particle types

    void DumpCouples() const;
      // Displays a list of registered couples

    const G4MCCIndexConversionTable* GetMCCIndexConversionTable() const;
      // Gives the pointer to the MCCIndexConversionTable

    const std::vector<G4double>* GetRangeCutsVector(std::size_t pcIdx) const;
    const std::vector<G4double>* GetEnergyCutsVector(std::size_t pcIdx) const;
  
    std::size_t GetTableSize() const;
      // Returns the size of the couple table

    const G4MaterialCutsCouple* GetMaterialCutsCouple(G4int i) const;
      // Returns the pointer to the couple

    const G4MaterialCutsCouple* GetMaterialCutsCouple(const G4Material* aMat,
                                            const G4ProductionCuts* aCut) const;
     // Returns the pointer to the couple

    G4int GetCoupleIndex(const G4MaterialCutsCouple* aCouple) const;
    G4int GetCoupleIndex(const G4Material* aMat,
                         const G4ProductionCuts* aCut) const;
      // Return the index of the couple.
      // -1 is returned if index is not found

    G4bool IsModified() const;
      // Returns TRUE if at least one production cut value is modified
  
    void PhysicsTableUpdated();
      // Resets the status of IsModified(). This method must be exclusively
      // used by the RunManager when physics tables are built

    G4ProductionCuts* GetDefaultProductionCuts() const;
      // Returns the default production cuts
 
    G4double ConvertRangeToEnergy(const G4ParticleDefinition* particle,
                                  const G4Material* material, 
                                        G4double range);
      // Gives energy corresponding to range value.
      // -1 is returned if particle or material is not found

    void ResetConverters();
      // Resets all range to energy converters  

    G4bool StoreCutsTable(const G4String& directory, 
                          G4bool ascii = false);
      // Stores cuts and material information in files under the
      // the specified directory

    G4bool RetrieveCutsTable(const G4String& directory,
                             G4bool ascii = false);
      // Retrieve material cut couple information 
      // in files under the specified directory

    G4bool CheckForRetrieveCutsTable(const G4String& directory, 
                                     G4bool ascii = false);
      // Checks stored material and cut values are consistent
      // with the current detector setup

    G4double* GetRangeCutsDoubleVector(std::size_t pcIdx) const;
    G4double* GetEnergyCutsDoubleVector(std::size_t pcIdx) const;
      // Methods for backward compatibility

    void SetEnergyCutVector(const std::vector<G4double>& cutE, std::size_t idx);
      // User defined cut vectors (idx < 4) range cut should be defined
      // to avoid inconsistency in physics

    void  SetVerboseLevel(G4int value);
    G4int GetVerboseLevel() const;
      // Control flag for output message
      //  0: Silent
      //  1: Warning message
      //  2: More

  protected:

    G4ProductionCutsTable();

    virtual G4bool StoreMaterialInfo(const G4String& directory, 
                                     G4bool ascii = false);
      // Stores material information in files under the specified directory

    virtual G4bool CheckMaterialInfo(const G4String& directory, 
                                     G4bool ascii = false);
      // Checks stored material is consistent with the current detector setup

    virtual G4bool StoreMaterialCutsCoupleInfo(const G4String& directory, 
                                               G4bool ascii = false);
      // Stores materialCutsCouple information in files under the
      // specified directory

    virtual G4bool CheckMaterialCutsCoupleInfo(const G4String& directory, 
                                               G4bool ascii = false);
      // Checks stored materialCutsCouple is consistent with
      // the current detector setup

    virtual G4bool StoreCutsInfo(const G4String& directory, 
                                 G4bool ascii = false);
      // Stores cut values information in files under the specified directory

    virtual G4bool  RetrieveCutsInfo(const G4String& directory,
                                     G4bool ascii = false);
      // Retrieves cut values information in files under the
      // specified directory

  private:

    void ScanAndSetCouple(G4LogicalVolume* aLV,
                          G4MaterialCutsCouple* aCouple,
                          G4Region* aRegion);

    G4bool IsCoupleUsedInTheRegion(const G4MaterialCutsCouple* aCouple,
                                   const G4Region* aRegion) const;


  private:

    static G4ProductionCutsTable* fProductionCutsTable;

    std::vector<G4MaterialCutsCouple*> coupleTable;
    std::vector<std::vector<G4double>*> rangeCutTable;
    std::vector<std::vector<G4double>*> energyCutTable;

    std::vector<G4double>* userEnergyCuts[4] = {nullptr, nullptr, nullptr, nullptr};

    G4RegionStore* fG4RegionStore = nullptr;
    G4VRangeToEnergyConverter* converters[NumberOfG4CutIndex]; 

    G4ProductionCuts* defaultProductionCuts = nullptr;

    G4MCCIndexConversionTable mccConversionTable;

    // These two vectors are for backward compatibility
    G4double* rangeDoubleVector[NumberOfG4CutIndex];
    G4double* energyDoubleVector[NumberOfG4CutIndex];

    enum { FixedStringLengthForStore = 32 }; 

    G4ProductionCutsTableMessenger* fMessenger = nullptr;
    G4int verboseLevel = 1;
    G4bool firstUse = true;
};

// ------------------
// Inline methods
// ------------------

inline 
const std::vector<G4double>*
G4ProductionCutsTable::GetRangeCutsVector(std::size_t pcIdx) const
{ 
  return rangeCutTable[pcIdx]; 
}

inline 
const std::vector<G4double>*
G4ProductionCutsTable::GetEnergyCutsVector(std::size_t pcIdx) const
{
  return energyCutTable[pcIdx]; 
}

inline 
std::size_t G4ProductionCutsTable::GetTableSize() const
{
  return coupleTable.size(); 
}

inline 
const G4MaterialCutsCouple*
G4ProductionCutsTable::GetMaterialCutsCouple(G4int i) const
{ 
  return coupleTable[std::size_t(i)]; 
}

inline 
G4bool G4ProductionCutsTable::IsModified() const
{ 
  if(firstUse) return true;
  for(auto itr=coupleTable.cbegin(); itr!=coupleTable.cend(); ++itr)
  { 
    if((*itr)->IsRecalcNeeded())
    {
      return true; 
    }
  }
  return false;
}
  
inline 
void G4ProductionCutsTable::PhysicsTableUpdated()
{ 
  for(auto itr=coupleTable.cbegin(); itr!=coupleTable.cend(); ++itr)
  { 
    (*itr)->PhysicsTableUpdated(); 
  }
}

inline
G4double*
G4ProductionCutsTable::GetRangeCutsDoubleVector(std::size_t pcIdx) const
{
  return rangeDoubleVector[pcIdx];
}

inline
G4double*
G4ProductionCutsTable::GetEnergyCutsDoubleVector(std::size_t pcIdx) const
{
  return energyDoubleVector[pcIdx];
}

inline
G4ProductionCuts* G4ProductionCutsTable::GetDefaultProductionCuts() const
{
  return defaultProductionCuts;
}

inline
G4bool G4ProductionCutsTable::IsCoupleUsedInTheRegion(
                                 const G4MaterialCutsCouple* aCouple,
                                 const G4Region* aRegion) const
{
  G4ProductionCuts* fProductionCut = aRegion->GetProductionCuts();
  auto mItr = aRegion->GetMaterialIterator();
  std::size_t nMaterial = aRegion->GetNumberOfMaterials();
  for(std::size_t iMate=0;iMate<nMaterial; ++iMate, ++mItr)
  {
    if(aCouple->GetMaterial()==(*mItr) &&
       aCouple->GetProductionCuts()==fProductionCut)
    {
      return true;
    }
  }
  return false;
}

inline
const G4MaterialCutsCouple* 
G4ProductionCutsTable::GetMaterialCutsCouple(const G4Material* aMat,
                                             const G4ProductionCuts* aCut) const
{
  for(auto cItr=coupleTable.cbegin(); cItr!=coupleTable.cend(); ++cItr)
  { 
    if((*cItr)->GetMaterial()!=aMat) continue;
    if((*cItr)->GetProductionCuts()==aCut) return (*cItr);
  }
  return nullptr;
}

inline
G4int
G4ProductionCutsTable::GetCoupleIndex(const G4MaterialCutsCouple* aCouple) const
{ 
  G4int idx = 0;
  for(auto cItr=coupleTable.cbegin(); cItr!=coupleTable.cend(); ++cItr)
  {
    if((*cItr)==aCouple) return idx;
    ++idx;
  }
  return -1;
}

inline
G4int G4ProductionCutsTable::GetCoupleIndex(const G4Material* aMat,
                                            const G4ProductionCuts* aCut) const
{
  const G4MaterialCutsCouple* aCouple = GetMaterialCutsCouple(aMat,aCut);
  return GetCoupleIndex(aCouple);
}

inline 
G4int G4ProductionCutsTable::GetVerboseLevel() const
{
  return verboseLevel;
}

inline
const G4MCCIndexConversionTable* 
G4ProductionCutsTable::GetMCCIndexConversionTable() const
{
  return &mccConversionTable;
}

#endif
