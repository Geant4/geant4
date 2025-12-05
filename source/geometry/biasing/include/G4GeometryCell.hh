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
// G4GeometryCell
//
// Class description:
//
// This class is used by scoring and importance sampling.
// It serves to address a "cell". A "cell" is somewhat related to a
// touchable. It is identified by a reference to a G4VPhysicalVolume
// and a number (replica number). Only simple replicas are supported.

// Author: Michael Dressel (CERN), 2002
// ----------------------------------------------------------------------
#ifndef G4GEOMETRYCELL_HH
#define G4GEOMETRYCELL_HH

#include "globals.hh"

class G4VPhysicalVolume;

/**
 * @brief G4GeometryCell is used for scoring and importance sampling.
 * It defines a "cell", which, similar to the concept of "touchable",
 + is identified by a reference to a G4VPhysicalVolume and a number
 * (replica number).
 */

class G4GeometryCell
{
  public:

    /**
     * Constructor. Initialises volume and replica number.
     *  @param[in] aVolume The name of the physical volume.
     *  @param[in] RepNum The associated replica number.
     */
    G4GeometryCell(const G4VPhysicalVolume& aVolume, G4int RepNum);

    /**
     * Copy constructor and assignment operator.
     */
    G4GeometryCell(const G4GeometryCell& rhs) = default;
    G4GeometryCell& operator=(const G4GeometryCell& rhs);

    /**
     * Default Destructor.
     */
    ~G4GeometryCell() = default;
  
    /**
     * Returns the physical volume of the cell.
     */
    const G4VPhysicalVolume& GetPhysicalVolume() const;

    /**
     * Returns the replica number of the cell.
     */
    G4int GetReplicaNumber() const;
 
  private:

    /** Pointer to the physical volume of the "cell"; treated as identifier. */
    const G4VPhysicalVolume* fVPhysicalVolume = nullptr;

    /** Replica number of the "cell". */
    G4int fRepNum = 0;
};

G4bool operator==(const G4GeometryCell& k1, const G4GeometryCell& k2);
G4bool operator!=(const G4GeometryCell& k1, const G4GeometryCell& k2);

#endif
