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
// $Id: G4IStore.hh 102994 2017-03-07 16:31:28Z gcosmo $
//
// ----------------------------------------------------------------------
// Class G4IStore
//
// Class description:
//
// An implementation of a "importance store" with the interface
// G4VIStore. See description in G4VIStore.
// This implementation uses G4GeometryCellImportance as the container
// to store the "cells" together with the importance values.
// Giving a cell the importance 0 is allowed as
// a flagging that no biasing should happen between this
// cell and it's neighbors
// If a cell is not known by the importance store no biasing
// should be applied between this cell and it's nighbors.

// Author: Michael Dressel (Michael.Dressel@cern.ch)
//
// Alex Howard (alexander.howard@cern.ch):
// Changed class to a `singleton', with access via the static method
// G4IStore::GetInstance().
//
// Member data:
//
//   static G4IStore* fInstance
//     - Ptr to the unique instance of class
// ----------------------------------------------------------------------
#ifndef G4IStore_hh
#define G4IStore_hh G4IStore_hh 

#include "G4VIStore.hh"
#include "G4GeometryCellImportance.hh"
#include "G4TransportationManager.hh"

class G4IStore : public G4VIStore
{

public:

  //  static G4ThreadLocal G4IStore* GetInstance();
  static G4IStore* GetInstance();
  // Return ptr to singleton instance of the class.
  //  static G4ThreadLocal G4IStore* GetInstance(G4String ParallelWorldName);
  static G4IStore* GetInstance(const G4String& ParallelWorldName);
  // Return ptr to singleton instance of the class.

protected:

  explicit G4IStore();
    // initialise the importance store for the given geometry
  explicit G4IStore(const G4String& ParallelWorldName);
    // initialise the importance store for the given geometry

  ~G4IStore();
    // destruct

public:  // with description

  virtual G4double GetImportance(const G4GeometryCell &gCell) const;
    // derive an importance value of a "cell" addressed by a G4GeometryCell
    // from the store.

  virtual G4bool IsKnown(const G4GeometryCell &gCell) const;
    // returns true if the gCell is in the store, else false 

  void Clear();

  void SetWorldVolume();
    // set a reference to the world volume of the 
    // "importance" geometry

  void SetParallelWorldVolume(G4String paraName);
    // set a reference to the parallel world volume of the 
    // "importance" geometry

  virtual const G4VPhysicalVolume& GetWorldVolume() const;
    // return a reference to the world volume of the 
    // "importance" geometry

  // virtual const G4VPhysicalVolume& GetParallelWorldVolume() const;
  //   // return a reference to the world volume of the 
  //   // "importance" geometry
  virtual const G4VPhysicalVolume* GetParallelWorldVolumePointer() const;
    // return a pointer to the world volume of the 
    // "importance" geometry

  void AddImportanceGeometryCell(G4double importance,
			   const G4GeometryCell &gCell);
  void AddImportanceGeometryCell(G4double importance,
			   const G4VPhysicalVolume &,
			   G4int aRepNum = 0);
    // Add a "cell" together with a importance value to the store.

  void ChangeImportance(G4double importance,
			const G4GeometryCell &gCell);
  void ChangeImportance(G4double importance,
			const G4VPhysicalVolume &,
			G4int aRepNum = 0);
    // Change a importance value of a "cell".

  G4double GetImportance(const G4VPhysicalVolume &,
			 G4int aRepNum = 0) const ;
  
private:

  G4bool IsInWorld(const G4VPhysicalVolume &) const;
  void SetInternalIterator(const G4GeometryCell &gCell) const;
  void Error(const G4String &m) const;
  
private:
 
  const G4VPhysicalVolume* fWorldVolume;
  //  const G4VPhysicalVolume* fParallelWorldVolume;
  //  G4bool fParaFlag;
  G4GeometryCellImportance fGeometryCelli;

  mutable G4GeometryCellImportance::const_iterator fCurrentIterator;

  static G4ThreadLocal G4IStore* fInstance;

#ifdef G4MULTITHREADED
  static G4Mutex IStoreMutex;
#endif

};

#endif
