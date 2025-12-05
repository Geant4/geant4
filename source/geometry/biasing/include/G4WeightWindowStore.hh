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
// G4WeightWindowStore
//
// Class description:
//
// Implementation of a weight window store according to the
// G4VWeightWindowStore interface.
// See also G4VWeightWindowStore.

// Author: Michael Dressel (CERN), 2003 - Created
//         Alex Howard (CERN), 2013 - Changed class to a 'singleton'
// ----------------------------------------------------------------------
#ifndef G4WEIGHTWINDOWSTORE_HH
#define G4WEIGHTWINDOWSTORE_HH

#include "G4VWeightWindowStore.hh"
#include "G4GeometryCellWeight.hh"
#include <set>
#include <vector>

/**
 * @brief G4WeightWindowStore is an concrete implementation of a weight window
 * store according to the G4VWeightWindowStore interface.
 */

class G4WeightWindowStore: public G4VWeightWindowStore
{
  public:

    /**
     * Returns a pointer to the singleton instance of the class.
     */
    static G4WeightWindowStore* GetInstance();

    /**
     * Returns a pointer to the singleton instance of the class, given
     * the name of the parallel world of reference.
     */
    static G4WeightWindowStore* GetInstance(const G4String& ParallelWorldName);

    /**
     * Derives a lower weight bound value of a "cell" addressed by a 
     * G4GeometryCell and the corresponding energy from the store.
     */
    G4double GetLowerWeight(const G4GeometryCell& gCell, 
                                  G4double partEnergy) const override;

    /**
     * Returns true if 'gCell' is in the store, else false.
     *  @param[in] gCell The cell of reference.
     *  @returns true if present in the store, false otherwise.
     */
    G4bool IsKnown(const G4GeometryCell &gCell) const override;

    /**
     * Clears the cells weights map.
     */
    void Clear();

    /**
     * Sets a reference to world volume of the "weightwindow" geometry.
     */
    void SetWorldVolume();

    /**
     * Sets a reference to parallel world volume of the "weightwindow" geometry.
     */
    void SetParallelWorldVolume(const G4String& paraName);

    /**
     * Returns a reference to the world volume of the "weightwindow" geometry.
     */
    const G4VPhysicalVolume& GetWorldVolume() const override;

    /**
     * Returns a pointer to the world volume of the "weightwindow" geometry.
     */
    const G4VPhysicalVolume* GetParallelWorldVolumePointer() const;

    /**
     * Adds lower weights. Only if general upper energy bounds have been set.
     */
    void AddLowerWeights(const G4GeometryCell& gCell,
                         const std::vector<G4double>& lowerWeights);

    /**
     * Sets upper energy - lower weight pairs for a cell.
     */
    void AddUpperEboundLowerWeightPairs(const G4GeometryCell& gCell,
                                  const G4UpperEnergyToLowerWeightMap& enWeMap);

    /**
     * Sets the energy bounds.
     */
    void SetGeneralUpperEnergyBounds(const std::set<G4double,
                                           std::less<G4double> >& enBounds);
  private:

    /**
     * Constructors. Initialise the weight window store for the given geometry.
     */
    explicit G4WeightWindowStore();
    explicit G4WeightWindowStore(const G4String& ParallelWorldName);

    /**
     * Default Destructor.
     */
    ~G4WeightWindowStore() override = default;

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

    std::set<G4double, std::less<G4double> > fGeneralUpperEnergyBounds;
    G4GeometryCellWeight fCellToUpEnBoundLoWePairsMap;
    mutable G4GeometryCellWeight::const_iterator fCurrentIterator;

    static G4ThreadLocal G4WeightWindowStore* fInstance;
};

#endif
