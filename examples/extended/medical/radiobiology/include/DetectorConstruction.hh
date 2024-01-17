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
/// \file radiobiology/include/DetectorConstruction.hh
/// \brief Definition of the RadioBio::DetectorConstruction class

#ifndef RadiobiologyDetectorConstruction_h
#define RadiobiologyDetectorConstruction_h 1

#include "G4SystemOfUnits.hh"
#include "G4ThreeVector.hh"
#include "G4VUserDetectorConstruction.hh"
#include "globals.hh"

class G4LogicalVolume;
class G4Material;

namespace RadioBio
{

// Forward declariation of other radiobiology classes
class DetectorMessenger;

/**
 * @brief Mandatory class for the construction of geometry.
 *
 * This class is used to create the world volume and keep track of some
 * dimensions. World volume is the only one created and saved here,
 * whereas voxels and sensitiveness are accounted for by another
 * class.
 *
 * @note Different class will use some information to this class such
 * as geometrical constraints on the simulation world volume.
 *
 */

class DetectorConstruction : public G4VUserDetectorConstruction
{
  public:
    DetectorConstruction();
    ~DetectorConstruction();

  public:
    /** @brief Method overriding pure virtual
     * G4VUserDetectorConstruction one */
    G4VPhysicalVolume* Construct() override;
    /** @brief Method overriding pure virtual
     *  G4VUserDetectorConstruction one */
    void ConstructSDandField() override;

    // Set cubical world size
    /** @brief Method to set the world volume
     *  to be cubical with given side */
    void SetSize(G4double);

    // Set world size for non-cubic boxes
    /** @brief Method to set the world volume
     *  to be a box with given (potentially different) sides */
    void SetSize(G4ThreeVector);
    /** @brief Method to set the X width of the world volume
     *  @param x X width of the box */
    void SetSizeX(G4double x);
    /** @brief Method to set the Y width of the world volume
     *  @param y Y width of the box */
    void SetSizeY(G4double y);
    /** @brief Method to set the Z width of the world volume
     *  @param z Z width of the box */
    void SetSizeZ(G4double z);

    /** @brief Method to set the world material
     *
     * @note there will not be any different material throughout the
     * whole simulation
     *
     * @param mat string containing the name of the material
     * according to NIST database */
    void SetMaterial(G4String mat);

  public:
    /** @brief Returns a pointer to the world physical volume */
    const G4VPhysicalVolume* GetWorld() { return fPBox; }

    G4double GetSizeX() const { return fBoxSizeX; }
    G4double GetSizeY() const { return fBoxSizeY; }
    G4double GetSizeZ() const { return fBoxSizeZ; }


    /** @brief Returns a pointer to the world material */
    G4Material* GetMaterial()
    {
      return fMaterial;
    }

    /** @brief Prints on screen some parameters for this class */
    void PrintParameters();

  private:
    // World physical and logical
    G4VPhysicalVolume* fPBox = nullptr;
    G4LogicalVolume* fLBox = nullptr;

    // World box sizes
    G4double fBoxSizeX = 1. * m;
    G4double fBoxSizeY = 1. * m;
    G4double fBoxSizeZ = 1. * m;

    // World box material
    G4Material* fMaterial = nullptr;

    // Pointer to the messenger
    DetectorMessenger* fDetectorMessenger = nullptr;

  private:
    /** @brief Private method to construct geometry */
    G4VPhysicalVolume* ConstructVolumes();
};

}  // namespace RadioBio

#endif  // DetectorConstruction_h
