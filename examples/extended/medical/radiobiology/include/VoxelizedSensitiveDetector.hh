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
/// \file radiobiology/include/VoxelizedSensitiveDetector.hh
/// \brief Definition of the RadioBio::VoxelizedSensitiveDetector class

#ifndef RadiobiologyVoxelizedSensitiveDetector_h
#define RadiobiologyVoxelizedSensitiveDetector_h 1

#include "G4ThreeVector.hh"

class G4Box;
class G4LogicalVolume;
class G4VPhysicalVolume;

namespace RadioBio
{

// Forward declariation of other radiobiology classes
class DetectorConstruction;
class VoxelizedSensitiveDetectorMessenger;

/**
 * @brief Singleton class performing the voxelization of the
 * world volume and tracking of a voxel index given the three-dimensional
 * coordinates
 *
 * This class is used to voxelize the world volume created in
 * DetectorConstruction class. It works only if the world already exists:
 * the constructor needs a pointer to DetectorConstruction.
 *
 * This class does not insert a further material, but uses the same
 * material saved in DetectorConstruction class. It is responsible
 * for the voxelization of the space and for the tracking of voxels
 * index.
 *
 * Internally, the voxel volume, mass and density are updated every
 * time something relevant changes in the geometry
 *
 * @note Although in this class the word \a "sensitive" appears,
 * this is class is not an implementation of the sensitive detector,
 * nevertheless it crates the environment for sensitiveness. It is
 * meant to be similar to a DetectorConstruction class, providing
 * voxelization of the world.
 */

class VoxelizedSensitiveDetector
{
    // This class has a 'singleton' structure

  private:
    /** @brief Private constructor using a pointer to
     * DetectorConstruction and the three dimensions for one voxel */
    VoxelizedSensitiveDetector(DetectorConstruction* det, double xWidth, double yWidth,
                               double zWidth);
    static VoxelizedSensitiveDetector* fInstance;

  public:
    /** @brief Static method to retrieve a pointer
     * to the only object existing in the simulation */
    static VoxelizedSensitiveDetector* GetInstance();

    /** @brief Static method to create the pointer
     * to the only object existing in the simulation
     *
     * @param det A pointer to the DetectorConstruction object
     * @param xWidth X dimension of a single voxel
     * @param yWidth Y dimension of a single voxel
     * @param zWidth Z dimension of a single voxel
     * */
    static VoxelizedSensitiveDetector* CreateInstance(DetectorConstruction* det, double xWidth,
                                                      double yWidth, double zWidth);

    /** Destructor */
    ~VoxelizedSensitiveDetector();

    // Setting of parameters
    /** @brief Method to set the voxel shape
     *  to be a box with given sides
     *
     * @param voxWidth ThreeVector containing the
     * three-dimensionals sides for the voxel
     */
    void SetVoxelWidth(G4ThreeVector voxWidth);

    /** @brief Method to set the voxel X width
     *
     * @param voxWidthX X width for the voxel
     */
    void SetVoxelWidthX(G4double voxWidthX);

    /** @brief Method to set the voxel Y width
     *
     * @param voxWidthY Y width for the voxel
     */
    void SetVoxelWidthY(G4double voxWidthY);

    /** @brief Method to set the voxel Z width
     *
     * @param voxWidthZ Z width for the voxel
     */
    void SetVoxelWidthZ(G4double voxWidthZ);

    /** @brief Method to update voxelized geometry.
     * If voxelized geometry is not built, it returns
     * immediately.
     * Otherwise, geometrical parameters are re-calculated and
     * voxelized geometry is destroyed, de-registered
     *  and built once again from scratch.
     *
     */
    void UpdateVoxelizedGeometry();

    // Initialize pointer to the world
    // Should be called only from DetectorConstruction
    /** @brief Method to properly initialize the pointer
     * to the World physical volume. It is meant
     * to be called **only** by DetectorConstruction
     * itself.
     *
     * @param pWorld Pointer to the world physical volume
     *
     */
    void InitializeWorldPtr(G4VPhysicalVolume* pWorld);

    // 'Construct' methods
    /** Method to construct all the volumes for the
     * voxelization of the geometry. */
    void Construct();

    /** Method to make the proper volume
     * sensitive to allow scoring. */
    void ConstructSD();

    // List of 'Get' methods
  public:
    /** Method to get the three vector
     * containing the dimensions for the voxel */
    G4ThreeVector GetVoxelWidth() const
    {
      return G4ThreeVector{fVoxelWidthX, fVoxelWidthY, fVoxelWidthZ};
    }

    /** Method to get X width of the voxel. */
    G4double GetVoxelWidthX() const
    {
      return fVoxelWidthX;
    }

    /** Method to get Y width of the voxel. */
    G4double GetVoxelWidthY() const
    {
      return fVoxelWidthY;
    }

    /** Method to get Z width of the voxel. */
    G4double GetVoxelWidthZ() const
    {
      return fVoxelWidthZ;
    }

    /** Method to get total voxel volume. */
    G4double GetVoxelVolume() const
    {
      return fVoxelVolume;
    }

    /** Method to get the voxel number along X axis. */
    G4double GetVoxelNumberAlongX() const
    {
      return fVoxelNumberAlongX;
    }

    /** Method to get the voxel number along Y axis. */
    G4double GetVoxelNumberAlongY() const
    {
      return fVoxelNumberAlongY;
    }

    /** Method to get the voxel number along Z axis. */
    G4double GetVoxelNumberAlongZ() const
    {
      return fVoxelNumberAlongZ;
    }

    /** Method to get the total voxel number. */
    G4int GetTotalVoxelNumber() const
    {
      return fTotalVoxelNumber;
    }

    /** Method to get the mass of a voxel.
     *
     * @note In this implementation, there is only
     * one material for the whole simulation. This
     * means that all the voxels (that already share the
     * same shape and dimensions) also share the same
     * mass and density
     */
    G4double GetVoxelMass() const
    {
      return fVoxelMass;
    }

    /** Method to get the density of a voxel.
     *
     * @note In this implementation, there is only
     * one material for the whole simulation. This
     * means that all the voxels (that already share the
     * same shape and dimensions) also share the same
     * mass and density
     */
    G4double GetVoxelDensity() const
    {
      return fVoxelDensity;
    }

    // To get the onedimensional voxel number given the three indexes
    // This function will be used from other classes as well, and
    // is defined only here for clarity purpose
    /** Method to get the absolute voxel index given its
     * indexes in the three dimensions.
     *
     * This helper function is heavily used throughout
     * the simulation and is here because the voxelized
     * detector is meant to know how many voxels there
     * are and how they are ordered.
     */
    G4int GetThisVoxelNumber(G4int x, G4int j, G4int k) const;

  private:
    /** Private method to construct the voxelized
     * detector
     */
    G4bool ConstructVoxelizedDetector();

    DetectorConstruction* fDetector = nullptr;

    G4double fVoxelWidthX = -1.;
    G4double fVoxelWidthY = -1.;
    G4double fVoxelWidthZ = -1.;

    G4double fVoxelVolume = -1.;

    G4int fVoxelNumberAlongX = 1;
    G4int fVoxelNumberAlongY = 1;
    G4int fVoxelNumberAlongZ = 1;

    G4int fTotalVoxelNumber = 1;
    G4double fVoxelMass = 0.;
    G4double fVoxelDensity = 0.;

    // Voxelized objects
    G4Box* fVoxelizedDetectorXDivision = nullptr;
    G4Box* fVoxelizedDetectorYDivision = nullptr;
    G4Box* fVoxelizedDetectorZDivision = nullptr;

    G4LogicalVolume* fVoxelizedDetectorXDivisionLog = nullptr;
    G4LogicalVolume* fVoxelizedDetectorYDivisionLog = nullptr;
    G4LogicalVolume* fVoxelizedDetectorZDivisionLog = nullptr;
    G4LogicalVolume* fSensitiveLogicalVolume = nullptr;

    G4VPhysicalVolume* fVoxelizedDetectorXDivisionPhys = nullptr;
    G4VPhysicalVolume* fVoxelizedDetectorYDivisionPhys = nullptr;
    G4VPhysicalVolume* fVoxelizedDetectorZDivisionPhys = nullptr;

    // Pointer to world
    G4LogicalVolume* fWorldLogical = nullptr;

    G4bool fIsBuilt = false;

    /** Private method to calculate
     * the total volume of a voxel given the
     * parameters saved inside the object
     */
    void UpdateVoxelVolume();

    /** Private method to calculate
     * the total voxel number given the
     * parameters saved inside the object
     */
    void CalculateVoxelNumber();

    /** Private method to slice
     * the world volume along the X
     * axis.
     */
    void ConstructXDivision();

    /** Private method to further slice
     * the X slices along the Y
     * axis, creating some "columns" of material.
     */
    void ConstructYDivision();

    /** Private method to further slice
     * the Y columns along the Z
     * axis, creating some boxes of material.
     *
     * This is the final step for voxelization
     * of space.
     */
    void ConstructZDivision();

    VoxelizedSensitiveDetectorMessenger* fVoxelizedSensitiveDetectorMessenger = nullptr;
};

}  // namespace RadioBio

#endif  // VoxelizedSensitiveDetector_h
