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
// G4VIStore
//
// Class description:
//
// An interface of an "importance store" used by importance sampling.
// It defines how a importance value together with a "cell" 
// (a G4VPhysicalVolume and a replica number) has to be added
// to the "importance store" and how a importance value can be derived 
// from the "importance store".

// Author: Michael Dressel (CERN), 2002
// ----------------------------------------------------------------------
#ifndef G4VISTORE_HH
#define G4VISTORE_HH

#include "globals.hh"

class G4GeometryCell;
class G4VPhysicalVolume;

/**
 * @brief G4VIStore is an interface of an "importance store" used by importance
 * sampling. It defines how an importance value together with a "cell" 
 * (a G4VPhysicalVolume and a replica number) has to be added to the
 * "importance store" and how a importance value can be derived from the
 * "importance store". 
 */

class  G4VIStore
{
  public:

    /**
     * Default Constructor and Destructor.
     */
    G4VIStore() = default;
    virtual ~G4VIStore() = default;

    /**
     * Returns the importance value of a "cell" from the store addressed
     * by 'gCell'.
     *  @param[in] gCell The cell of reference.
     *  @returns The associated importance weight.
     */
    virtual G4double GetImportance(const G4GeometryCell& gCell) const = 0;

    /**
     * Returns true if 'gCell' is in the store, else false.
     *  @param[in] gCell The cell of reference.
     *  @returns true if present in the store, false otherwise.
     */
    virtual G4bool IsKnown(const G4GeometryCell& gCell) const = 0;

    /**
     * Returns a reference to the world volume of the "importance" geometry.
     */
    virtual const G4VPhysicalVolume& GetWorldVolume() const = 0;
};

#endif
