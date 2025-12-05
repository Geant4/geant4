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

// Author: Michael Dressel (CERN), 2002 - Created
//         Alex Howard (CERN), 2013 - Changed class to a 'singleton'
// ----------------------------------------------------------------------
#ifndef G4ISTORE_HH
#define G4ISTORE_HH

#include "G4VIStore.hh"
#include "G4GeometryCellImportance.hh"
#include "G4TransportationManager.hh"

/**
 * @brief G4IStore is a concrete implementation of an "importance store", as
 * derived from G4VIStore. It is a singleton, using G4GeometryCellImportance
 * as the container to store the "cells" together with the importance values.
 * Giving a cell, the importance 0 is allowed as a flagging that no biasing
 * should happen between this cell and its neighbors.
 * If a cell is not known by the importance store no biasing should be
 * applied between this cell and its nighbors.
 */

class G4IStore : public G4VIStore
{
  public:

    /**
     * Returns a pointer to the singleton instance of the class.
     */
    static G4IStore* GetInstance();

    /**
     * Returns a pointer to the singleton instance of the class, given
     * the name of the parallel world of reference.
     */
    static G4IStore* GetInstance(const G4String& ParallelWorldName);

    /**
     * Returns the importance value of a "cell" from the store addressed
     * by 'gCell'.
     *  @param[in] gCell The cell of reference.
     *  @returns The associated importance weight.
     */
    G4double GetImportance(const G4GeometryCell& gCell) const override;

    /**
     * Returns true if 'gCell' is in the store, else false.
     *  @param[in] gCell The cell of reference.
     *  @returns true if present in the store, false otherwise.
     */
    G4bool IsKnown(const G4GeometryCell& gCell) const override;

    /**
     * Clears the cells importance store.
     */
    void Clear();

    /**
     * Sets a reference to world volume of the "importance" geometry.
     */
    void SetWorldVolume();

    /**
     * Sets a reference to parallel world volume of the "importance" geometry.
     */
    void SetParallelWorldVolume(const G4String& paraName);

    /**
     * Returns a reference to the world volume of the "importance" geometry.
     */
    const G4VPhysicalVolume& GetWorldVolume() const override;

    /**
     * Returns a pointer to the world volume of the "importance" geometry.
     */
    const G4VPhysicalVolume* GetParallelWorldVolumePointer() const;

    /**
     * Methods to add a "cell" together with an importance value to the store.
     */
    void AddImportanceGeometryCell(G4double importance,
                             const G4GeometryCell &gCell);
    void AddImportanceGeometryCell(G4double importance,
                             const G4VPhysicalVolume &,
                             G4int aRepNum = 0);

    /**
     * Methods to change an importance value of a "cell".
     */
    void ChangeImportance(G4double importance, const G4GeometryCell& gCell);
    void ChangeImportance(G4double importance, const G4VPhysicalVolume&,
                          G4int aRepNum = 0);

    /**
     * Returns the importance weight, given the volume and replica number.
     */
    G4double GetImportance(const G4VPhysicalVolume& vol, G4int rpNum = 0) const;
  
  private:

    /**
     * Constructors. Initialising the importance store for the given geometry.
     */
    explicit G4IStore();
    explicit G4IStore(const G4String& ParallelWorldName);

    /**
     * Default Destructor.
     */
    ~G4IStore() override = default;

    /**
     * Internal utilities.
     */
    G4bool IsInWorld(const G4VPhysicalVolume&) const;
    void SetInternalIterator(const G4GeometryCell& gCell) const;

    /**
     * Internal logger.
     */
    void Error(const G4String& m) const;

  private:
 
    const G4VPhysicalVolume* fWorldVolume = nullptr;
    G4GeometryCellImportance fGeometryCelli;

    mutable G4GeometryCellImportance::const_iterator fCurrentIterator;

    static G4ThreadLocal G4IStore* fInstance;
};

#endif
