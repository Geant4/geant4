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
// G4IStore
//
// Class description:
//
// An implementation of an "importance store" with the interface
// G4VIStore. See description in G4VIStore.
// This implementation uses G4GeometryCellImportance as the container
// to store the "cells" together with the importance values.
// Giving a cell the importance 0 is allowed as a flagging that no biasing
// should happen between this cell and its neighbors.
// If a cell is not known by the importance store no biasing should be
// applied between this cell and its nighbors.

// Author: Michael Dressel (CERN), 2002
// Modified: Alex Howard (CERN), 2013 - Changed class to a 'singleton'
// ----------------------------------------------------------------------
#ifndef G4ISTORE_HH
#define G4ISTORE_HH 1

#include "G4VIStore.hh"
#include "G4GeometryCellImportance.hh"
#include "G4TransportationManager.hh"

class G4IStore : public G4VIStore
{
  public:

    static G4IStore* GetInstance();
      // return ptr to singleton instance of the class.

    static G4IStore* GetInstance(const G4String& ParallelWorldName);
      // return ptr to singleton instance of the class.

    virtual G4double GetImportance(const G4GeometryCell& gCell) const;
      // derive an importance value of a "cell" addressed by a
      // G4GeometryCell from the store.

    virtual G4bool IsKnown(const G4GeometryCell& gCell) const;
      // returns true if the gCell is in the store, else false 

    void Clear();

    void SetWorldVolume();
      // set a reference to the world volume of the "importance" geometry

    void SetParallelWorldVolume(const G4String& paraName);
      // set a reference to parallel world volume of the "importance" geometry

    virtual const G4VPhysicalVolume& GetWorldVolume() const;
      // return a reference to the world volume of the "importance" geometry

    virtual const G4VPhysicalVolume* GetParallelWorldVolumePointer() const;
      // return a pointer to the world volume of the "importance" geometry

    void AddImportanceGeometryCell(G4double importance,
                             const G4GeometryCell &gCell);
    void AddImportanceGeometryCell(G4double importance,
                             const G4VPhysicalVolume &,
                             G4int aRepNum = 0);
      // add a "cell" together with a importance value to the store.

    void ChangeImportance(G4double importance, const G4GeometryCell& gCell);
    void ChangeImportance(G4double importance, const G4VPhysicalVolume&,
                          G4int aRepNum = 0);
      // change a importance value of a "cell".

    G4double GetImportance(const G4VPhysicalVolume&, G4int aRepNum = 0) const;
  
  protected:

    explicit G4IStore();
      // initialise the importance store for the given geometry
    explicit G4IStore(const G4String& ParallelWorldName);
      // initialise the importance store for the given geometry

    ~G4IStore();
      // destructor

  private:

    G4bool IsInWorld(const G4VPhysicalVolume&) const;
    void SetInternalIterator(const G4GeometryCell& gCell) const;
    void Error(const G4String& m) const;

  private:
 
    const G4VPhysicalVolume* fWorldVolume = nullptr;
    G4GeometryCellImportance fGeometryCelli;

    mutable G4GeometryCellImportance::const_iterator fCurrentIterator;

    static G4ThreadLocal G4IStore* fInstance;

#ifdef G4MULTITHREADED
    static G4Mutex IStoreMutex;
#endif

};

#endif
