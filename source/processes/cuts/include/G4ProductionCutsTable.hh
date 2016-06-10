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
// $Id: G4ProductionCutsTable.hh 70369 2013-05-29 14:59:24Z gcosmo $
//
// 
// ------------------------------------------------------------
//	GEANT 4 class header file 
//
// Class Description
//  G4ProductionCutsTable is a static singleton class of a table of
//  G4ProductionCuts objects. This class also manages tables of
//  production cut and energy cut for each particle type.
//
// ------------------------------------------------------------
//   First Implementation          05 Oct. 2002  M.Asai    
//
//   Modified                      03 Feb 2004 H.Kurashige
//    Modify RetrieveCutsTable to allow ordering of materials and 
//    couples can be different from one in file (i.e. at storing)
//   Modified                      20 Aug. 2004 H.Kurashige
//    Modify RetrieveCutsTable to allow materials and 
//    couples can be different from one in file (i.e. at storing)
//   Modified                      2 Mar. 2008 H.Kurashige
//    add messenger
// ------------------------------------------------------------

#ifndef G4ProductionCutsTable_h 
#define G4ProductionCutsTable_h 1

class G4RegionStore;
class G4VRangeToEnergyConverter;
class G4LogicalVolume;
class G4VPhysicalVolume;
class G4ProductionCuts;

class G4ProductionCutsTableMessenger;

#include "globals.hh"
#include <cmath>
#include "G4ios.hh"
#include <vector>
#include "G4MaterialCutsCouple.hh"
#include "G4MCCIndexConversionTable.hh"
#include "G4Region.hh"


class G4ProductionCutsTable  
{
  public: // with description
    static G4ProductionCutsTable* GetProductionCutsTable();
    // This static method returns the singleton pointer of this class object.
    // At the first invokation of this method, the singleton object is instantiated.

  protected:
    G4ProductionCutsTable();
  private:
    G4ProductionCutsTable(const G4ProductionCutsTable& right);

  public:
    virtual ~G4ProductionCutsTable();

  public: // with description
    void UpdateCoupleTable(G4VPhysicalVolume* currentWorld);
    // This method triggers an update of the table of G4ProductionCuts objects.

    void SetEnergyRange(G4double lowedge, G4double highedge);
    // This method sets the limits of energy cuts for all particles.

    G4double GetLowEdgeEnergy() const;
    G4double GetHighEdgeEnergy() const;
    // These methods get the limits of energy cuts for all particles.

    //  get/set max cut energy of RangeToEnergy Converter for all particle type
    G4double GetMaxEnergyCut();
    void SetMaxEnergyCut(G4double value);


    void DumpCouples() const;
    // Display a list of registored couples

    const G4MCCIndexConversionTable* GetMCCIndexConversionTable() const;
    // gives the pointer to the MCCIndexConversionTable

  private:

   static G4ProductionCutsTable* fG4ProductionCutsTable;

   typedef std::vector<G4MaterialCutsCouple*> G4CoupleTable;
   typedef std::vector<G4MaterialCutsCouple*>::const_iterator CoupleTableIterator;
   typedef std::vector<G4double> G4CutVectorForAParticle;
   typedef std::vector<G4CutVectorForAParticle*> G4CutTable;
   G4CoupleTable coupleTable;
   G4CutTable rangeCutTable;
   G4CutTable energyCutTable;

   G4RegionStore* fG4RegionStore;
   G4VRangeToEnergyConverter* converters[NumberOfG4CutIndex]; 

   G4ProductionCuts* defaultProductionCuts;

   G4MCCIndexConversionTable mccConversionTable;

// These two vectors are for the backward comparibility
   G4double* rangeDoubleVector[NumberOfG4CutIndex];
   G4double* energyDoubleVector[NumberOfG4CutIndex];

  public: 
   const std::vector<G4double>* GetRangeCutsVector(size_t pcIdx) const;
   const std::vector<G4double>* GetEnergyCutsVector(size_t pcIdx) const;

// These two vectors are for the backward comparibility
   G4double* GetRangeCutsDoubleVector(size_t pcIdx) const;
   G4double* GetEnergyCutsDoubleVector(size_t pcIdx) const;
  
  public: // with description
   size_t GetTableSize() const;
     // This method returns the size of the couple table.

   const G4MaterialCutsCouple* GetMaterialCutsCouple(G4int i) const;
    // This method returns the pointer to the couple.

   const G4MaterialCutsCouple*
     GetMaterialCutsCouple(const G4Material* aMat,
                           const G4ProductionCuts* aCut) const;
    // This method returns the pointer to the couple.

   G4int GetCoupleIndex(const G4MaterialCutsCouple* aCouple) const;
   G4int GetCoupleIndex(const G4Material* aMat,
                           const G4ProductionCuts* aCut) const;
    // These methods return the index of the couple.
    // -1 is returned if index is not found.

   G4bool IsModified() const;
    // This method returns TRUE if at least one production cut value is modified.
  
   void PhysicsTableUpdated();
    // This method resets the status of IsModified(). This method must
    // be exclusively used by RunManager when physics tables are built.

   G4ProductionCuts* GetDefaultProductionCuts() const;
    // This method returns the default production cuts.
 
   G4double ConvertRangeToEnergy(const G4ParticleDefinition* particle,
				 const G4Material*           material, 
				 G4double                    range    );
    // This method gives energy corresponding to range value  
    // 
    // -1 is returned if particle or material is not found.

  void ResetConverters();
    // reset all Range To Energy Converters  
    
  private:
    void ScanAndSetCouple(G4LogicalVolume* aLV,
			  G4MaterialCutsCouple* aCouple,
			  G4Region* aRegion);

    bool IsCoupleUsedInTheRegion(const G4MaterialCutsCouple* aCouple,
				 const G4Region* aRegion) const;

  public: // with description
   // Store cuts and material information in files under the specified directory.
  G4bool  StoreCutsTable(const G4String& directory, 
			 G4bool          ascii = false);
  
  // Retrieve material cut couple information 
  //  in files under the specified directory.
  G4bool  RetrieveCutsTable(const G4String& directory,
			    G4bool          ascii = false);
  
  // check stored material and cut values are consistent with the current detector setup. 
  G4bool CheckForRetrieveCutsTable(const G4String& directory, 
				   G4bool          ascii = false);

  protected:

  // Store material information in files under the specified directory.
  virtual G4bool  StoreMaterialInfo(const G4String& directory, 
				    G4bool          ascii = false);

  // check stored material is consistent with the current detector setup. 
  virtual G4bool  CheckMaterialInfo(const G4String& directory, 
				    G4bool          ascii = false);
  
  // Store materialCutsCouple information in files under the specified directory.
  virtual G4bool  StoreMaterialCutsCoupleInfo(const G4String& directory, 
				    G4bool          ascii = false);

  // check stored materialCutsCouple is consistent with the current detector setup. 
  virtual G4bool  CheckMaterialCutsCoupleInfo(const G4String& directory, 
				    G4bool          ascii = false);

  // Store cut values information in files under the specified directory.
  virtual G4bool  StoreCutsInfo(const G4String& directory, 
                                G4bool          ascii = false);
  
 // Retrieve cut values information in files under the specified directory.
  virtual G4bool  RetrieveCutsInfo(const G4String& directory,
                                   G4bool          ascii = false);

  private:
   G4bool firstUse;
   enum { FixedStringLengthForStore = 32 }; 

  public: // with description  
      void  SetVerboseLevel(G4int value);
      G4int GetVerboseLevel() const;
      // controle flag for output message
      //  0: Silent
      //  1: Warning message
      //  2: More

  private:
    G4int verboseLevel;
    G4ProductionCutsTableMessenger* fMessenger;

};

inline 
 const std::vector<G4double>* G4ProductionCutsTable::GetRangeCutsVector(size_t pcIdx) const
{ 
  return rangeCutTable[pcIdx]; 
}

inline 
 const std::vector<G4double>* G4ProductionCutsTable::GetEnergyCutsVector(size_t pcIdx) const
{
 return energyCutTable[pcIdx]; 
}

inline 
 size_t G4ProductionCutsTable::GetTableSize() const
{
  return coupleTable.size(); 
}

inline 
 const G4MaterialCutsCouple* G4ProductionCutsTable::GetMaterialCutsCouple(G4int i) const
{ 
  return coupleTable[size_t(i)]; 
}

inline 
 G4bool G4ProductionCutsTable::IsModified() const
{ 
  if(firstUse) return true;
  for(G4ProductionCutsTable::CoupleTableIterator itr=coupleTable.begin();
    itr!=coupleTable.end();itr++){ 
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
  for(G4ProductionCutsTable::CoupleTableIterator itr=coupleTable.begin();itr!=coupleTable.end();itr++){ 
    (*itr)->PhysicsTableUpdated(); 
  }
}

inline
 G4double* G4ProductionCutsTable::GetRangeCutsDoubleVector(size_t pcIdx) const
{ return rangeDoubleVector[pcIdx]; }

inline
 G4double* G4ProductionCutsTable::GetEnergyCutsDoubleVector(size_t pcIdx) const
{ return energyDoubleVector[pcIdx]; }

inline
 G4ProductionCuts* G4ProductionCutsTable::GetDefaultProductionCuts() const
{ return defaultProductionCuts; }

inline
bool G4ProductionCutsTable::IsCoupleUsedInTheRegion(
                                 const G4MaterialCutsCouple* aCouple,
                                 const G4Region* aRegion) const
{
  G4ProductionCuts* fProductionCut = aRegion->GetProductionCuts();
  std::vector<G4Material*>::const_iterator mItr = aRegion->GetMaterialIterator();
  size_t nMaterial = aRegion->GetNumberOfMaterials();
  for(size_t iMate=0;iMate<nMaterial;iMate++, mItr++){
    if(aCouple->GetMaterial()==(*mItr) &&
       aCouple->GetProductionCuts()==fProductionCut){
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
  for(CoupleTableIterator cItr=coupleTable.begin();cItr!=coupleTable.end();cItr++)
  { 
    if((*cItr)->GetMaterial()!=aMat) continue;
    if((*cItr)->GetProductionCuts()==aCut) return (*cItr);
  }
  return 0;
}

inline
G4int G4ProductionCutsTable::GetCoupleIndex(const G4MaterialCutsCouple* aCouple) const
{ 
  G4int idx = 0;
  for(CoupleTableIterator cItr=coupleTable.begin();cItr!=coupleTable.end();cItr++)
  {
    if((*cItr)==aCouple) return idx;
    idx++;
  }
  return -1;
}

inline
G4int G4ProductionCutsTable:: GetCoupleIndex(const G4Material* aMat,
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






