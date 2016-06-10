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
// $Id: G4WeightWindowStore.hh 88801 2015-03-10 14:38:27Z gcosmo $
//
// ----------------------------------------------------------------------
// Class G4WeightWindowStore
//
// Class description:
//
// Implementation of a weight window store according to the
// G4VWeightWindowStore interface.
// See G4VWeightWindowStore.

// Author: Michael Dressel (Michael.Dressel@cern.ch)
//
// Alex Howard (alexander.howard@cern.ch):
// Changed class to a `singleton', with access via the static method
// G4WeightWindowStore::GetInstance().
//
// Member data:
//
//   static G4WeightWindowStore* fInstance
//     - Ptr to the unique instance of class
//
// ----------------------------------------------------------------------
#ifndef G4WeightWindowStore_hh
#define G4WeightWindowStore_hh G4WeightWindowStore_hh 

#include "G4VWeightWindowStore.hh"
#include "G4GeometryCellWeight.hh"
#include <set>
#include <vector>

class G4WeightWindowStore: public G4VWeightWindowStore
{

public:  // with description

  static G4WeightWindowStore* GetInstance();
  // Return ptr to singleton instance of the class.

  static G4WeightWindowStore* GetInstance(const G4String& ParallelWorldName);
  // Return ptr to singleton instance of the class.

protected:

  explicit G4WeightWindowStore();
    // initialise the weight window store for the given geometry
  explicit G4WeightWindowStore(const G4String& ParallelWorldName);
    // initialise the weight window store for the given geometry

  ~G4WeightWindowStore();
    // destructor

public:  // with description

  virtual G4double GetLowerWeight(const G4GeometryCell &gCell, 
			                G4double partEnergy) const;
    // derive a lower weight bound value of a "cell" addresed by a 
    // G4GeometryCell and the coresponding energy from the store.

  virtual G4bool IsKnown(const G4GeometryCell &gCell) const;
    // returns true if the gCell is in the store, else false 

  void Clear();

  void SetWorldVolume();
    // set a pointer to the world volume of the 
    // weightwindow geometry
  void SetParallelWorldVolume(G4String paraName);
    // set a pointer to the parallel world volume of the 
    // weightwindow geometry


  virtual const G4VPhysicalVolume &GetWorldVolume() const;
    // return a reference to the world volume of the 
    // weightwindow geometry
  virtual const G4VPhysicalVolume* GetParallelWorldVolumePointer() const;
    // return a pointer to the parallel world volume of the 
    // weightwindow geometry

  void AddLowerWeights(const G4GeometryCell &gCell,
		       const std::vector<G4double> &lowerWeights);
    // add lower weights. Only if general upper energy bounds have 
    // been set
 
  void AddUpperEboundLowerWeightPairs(const G4GeometryCell &gCell,
				      const G4UpperEnergyToLowerWeightMap&
				      enWeMap);
    // set upper energy - lower weight pairs for a cell

  void SetGeneralUpperEnergyBounds(const std::set<G4double, std::less<G4double> > &enBounds);
  
private:

  G4bool IsInWorld(const G4VPhysicalVolume &) const;
  void Error(const G4String &m) const;
  void SetInternalIterator(const G4GeometryCell &gCell) const;
  
private:

  const G4VPhysicalVolume* fWorldVolume;  
  //  G4bool fParaFlag;

  std::set<G4double, std::less<G4double> > fGeneralUpperEnergyBounds;
  G4GeometryCellWeight fCellToUpEnBoundLoWePairsMap;
  mutable G4GeometryCellWeight::const_iterator fCurrentIterator;

  static G4ThreadLocal G4WeightWindowStore* fInstance;
  
};

#endif
