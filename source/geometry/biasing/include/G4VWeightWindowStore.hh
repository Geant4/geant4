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
// G4VWeightWindowStore
//
// Class description:
//
// Interface class for a weight window store. It defines how the lower 
// weight window bound can be obtained from a weight window store.

// Author: Michael Dressel (CERN), 2003
// ----------------------------------------------------------------------
#ifndef G4VWEIGHTWINDOWSTORE_HH
#define G4VWEIGHTWINDOWSTORE_HH

#include "globals.hh"

class G4GeometryCell;
class G4VPhysicalVolume;

/**
 * @brief G4VWeightWindowStore is an interface class for a weight window store.
 * It defines how the lower weight window bound can be obtained from a weight
 * window store.
 */

class G4VWeightWindowStore
{
  public:

    /**
     * Default Constructor and Destructor.
     */
    G4VWeightWindowStore() = default;
    virtual ~G4VWeightWindowStore() = default;

    /**
     * Computes a lower weight bound value of a "cell" addressed by a 
     * G4GeometryCell and the corresponding energy from the store.
     */
    virtual G4double GetLowerWeight(const G4GeometryCell& gCell, 
                                          G4double partEnergy) const = 0;

    /**
     * Returns true if the gCell is in the store, else false.
     */
    virtual G4bool IsKnown(const G4GeometryCell& gCell) const = 0;

    /**
     * Returns a reference to the wolrd volume of the geometry.
     */
    virtual const G4VPhysicalVolume& GetWorldVolume() const = 0;
};

#endif
