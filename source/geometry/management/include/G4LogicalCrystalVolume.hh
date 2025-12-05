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
// G4LogicalCrystalVolume
//
// Class description:
//
// Specialised logical volume for the description of crystals.
// The object can be created only with an extended material
// with a crystal extension. The class handles the orientation
// of the crystal axes wrt the solid to which the crystal is attached.

// Author: Enrico Bagli (Ferrara Univ.), 21.04.2016
// ---------------------------------------------------------------------------
#ifndef G4LOGICALCRYSTALVOLUME_HH
#define G4LOGICALCRYSTALVOLUME_HH

#include <vector>
#include "G4LogicalVolume.hh"

class G4CrystalExtension;
class G4ExtendedMaterial;

/**
 * @brief G4LogicalCrystalVolume is a specialised logical volume for the
 * description of crystals. The object can be created only with an extended
 * material with a crystal extension. The class handles the orientation
 * of the crystal axes wrt the solid to which the crystal is attached.
 */

class G4LogicalCrystalVolume : public G4LogicalVolume
{
  public:

    /**
     * Constructor for G4LogicalCrystalVolume.
     *  @param[in] pSolid Pointer to the associated solid primitive.
     *  @param[in] pMaterial Pointer to the associated extended crystal material.
     *  @param[in] name The volume name.
     *  @param[in] pFieldMgr Pointer to optional magnetic field manager.
     *  @param[in] pSDetector Pointer to optional associated sensitive detector.
     *  @param[in] pULimits Pointer to optional user limits.
     *  @param[in] optimise Flag to enable/disable optimisation structure.
     *  @param[in] h Miller orientation parameter h.
     *  @param[in] k Miller orientation parameter k.
     *  @param[in] l Miller orientation parameter l.
     *  @param[in] rot Miller rotation parameter.
     */
    G4LogicalCrystalVolume(G4VSolid* pSolid,
                           G4ExtendedMaterial* pMaterial,
                           const G4String& name,
                           G4FieldManager* pFieldMgr = nullptr,
                           G4VSensitiveDetector* pSDetector = nullptr,
                           G4UserLimits* pULimits = nullptr,
                           G4bool optimise = true,
                           G4int h = 0,
                           G4int k = 0,
                           G4int l = 0,
                           G4double rot = 0.0);

    /**
     * Destructor.
     */
    ~G4LogicalCrystalVolume() override;
    
    /**
     * Returns true as it is not a base-class object.
     */
    inline G4bool IsExtended() const override { return true; }

    /**
     * Sets the physical lattice orientation, relative to G4VSolid coordinates.
     * Miller orientation aligns lattice normal (hkl) with geometry +Z.
     */
    void SetMillerOrientation(G4int h, G4int k, G4int l, G4double rot = 0.0);
    
    /**
     * Methods to rotate input vector between lattice and solid orientations.
     *  @returns The new vector value for convenience.
     */
    const G4ThreeVector& RotateToLattice(G4ThreeVector& dir) const;
    const G4ThreeVector& RotateToSolid(G4ThreeVector& dir) const;

    /**
     * Returns a pointer to the crystal extension object.
     */
    const G4CrystalExtension* GetCrystal() const;

    /**
     * Calls through to get crystal basis vectors.
     */
    const G4ThreeVector& GetBasis(G4int i) const;
    
    /**
     * Setter for verbosity level.
     */
    inline void SetVerbose(G4int aInt) { verboseLevel = aInt; }

    /**
     * Returns true if the logical volume 'aLV' is a lattice.
     */
    static G4bool IsLattice(G4LogicalVolume* aLV);

  private:

    /** Rotates geometry into lattice frame. */
    G4RotationMatrix fOrient;
    G4RotationMatrix fInverse;

    /** Cached Miller indices for dump. */
    G4int hMiller = 1, kMiller = 1, lMiller = 0;
    G4double fRot = 0.0;
    
    /** Verbosity level. */
    G4int verboseLevel = 0;

    static std::vector<G4LogicalVolume*> fLCVvec;
};

#endif
